"""
HAT_target_comparison.py
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
    against both targets both raw and with the scoring target's LOESS
    treatment (raw means D1-10, 10-domain LOESS beyond), the form the runs are
    graded in.

TWO MODEL SETS, ONE SUBFOLDER EACH (Hannah has not chosen the target)
    ends_solved_on_coastsat/   the matrix edgeBE run, end domains solved
                               against the CoastSat target
    ends_solved_on_duneline/   the 09-18 dune edge-solve run (mean3), end
                               domains solved against the dune line
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

The loaders are HAT_rate_windows.py's (imported), so the observations and
runs are exactly the ones model_vs_observed draws as rates.

USAGE
    python scripts/analyze_output/compare_runs/HAT_target_comparison.py
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
import HAT_rate_windows as rw  # noqa: E402  (loaders, runs, style constants)

_REPO = rw._REPO
sys.path.insert(0, str(_REPO / "scripts" / "input_prep" / "5-scr" / "duneline_vs_coastsat"))
from duneline_vs_coastsat import beach_width_handles, shade_beach_width  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.legend_handler import HandlerTuple  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import MultipleLocator  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    COMPARISONS_ROOT, DOMAIN_AXIS_LABEL, INK, _title, apply_style, caption,
    figsize, save,
)

obs = rw.obs
ROOT_DIR = COMPARISONS_ROOT / "target_comparison"
# The CoastSat target (2026-09-19, Hannah): the FULL-PERIOD 1996-2024 LRR in
# both windows, with runs whose ends were solved against it
# (experiments/2026-09-19-edgesolve-lrr1996_2024). The sub-period version
# (each window's own LRR, as the runner grades) is kept for the record.
CS_MODES = {"full": "coastsat_full_period_lrr", "subperiod": "coastsat_subperiod_lrr"}
FULL_WINDOW = (1996, 2024)
FULL_SOLVE_DIR = rw.RAW_RUNS / "experiments" / "2026-09-19-edgesolve-lrr1996_2024"
CS_MODE = "full"
OUT_DIR = ROOT_DIR / CS_MODES[CS_MODE]


def cs_label():
    return ("CoastSat target (1996–2024 LRR × 14 yr)" if CS_MODE == "full"
            else "CoastSat target (LRR × 14 yr)")


def cs_clause():
    return ("the full-period 1996–2024 linear regression rate (the same rate in both "
            "windows)" if CS_MODE == "full" else
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
        d = next((FULL_SOLVE_DIR / tag.split("/", 1)[1] / "{}_{}".format(*w) / "edgeBE").glob("HAT_*"))
        runs[w] = (d.name, tag)
    return runs
WINDOWS = [(1996, 2010), (2010, 2024)]
MODEL_SETS = {"coastsat": "ends_solved_on_coastsat",
              rw.MAIN_DUNE: "ends_solved_on_duneline"}
MODEL_LABEL = {"coastsat": "ends solved on CoastSat",
               rw.MAIN_DUNE: "ends solved on the dune line"}
LW = 1.1
LW_MODEL = 1.4
Y_LABEL = "Net change in position (m)"


def window_values(o, mdfs):
    """Per domain, metres over the model period: both targets raw and LOESS,
    and each model set's net change."""
    s, e = o.window
    years = e - s
    df = pd.DataFrame({"domain_number": np.arange(1, rw.N + 1)})
    df["model_years"] = years
    df["dune_interval_yr"] = o.meta["interval_yr"]
    cs, cs_t = CS_SOURCE.get(o.window, (o.coastsat, o.coastsat_target))
    df["coastsat_lrr_m_yr"] = cs["mean_lrr"].to_numpy(float)
    df["coastsat_target_m"] = df["coastsat_lrr_m_yr"] * years
    df["coastsat_target_loess_m"] = cs_t["target_lrr_m_yr"].to_numpy(float) * years
    df["dune_rate_m_yr"] = o.endpoint["mean_lrr"].to_numpy(float)
    df["dune_measured_change_m"] = df["dune_rate_m_yr"] * o.meta["interval_yr"]
    df["dune_target_m"] = df["dune_rate_m_yr"] * years
    df["dune_target_loess_m"] = o.endpoint_target["target_lrr_m_yr"].to_numpy(float) * years
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
                            ("coastsat_loess", "coastsat_target_loess_m"),
                            ("duneline", "dune_target_m"),
                            ("duneline_loess", "dune_target_loess_m")):
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
    obs.draw_shoals(ax, label=label)
    fills = obs.fills_in(*o.window)
    if fills:
        obs.draw_fills(ax, fills, half)


def figure(observations, frames, key, folder, half, skill_df):
    fig, axes = plt.subplots(len(WINDOWS), 1, sharex=True, sharey=True,
                             constrained_layout=True, figsize=figsize("double", height=5.6))
    for i, (ax, o) in enumerate(zip(axes, observations)):
        draw(ax, o, frames[o.window], folder, half, label=(i == 0))
        ax.yaxis.set_major_locator(MultipleLocator(20.0 if half > 60 else 10.0))
        _title(ax, i, "{}–{}".format(*o.window))
    axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)
    fig.supylabel(Y_LABEL, fontsize=9)
    handles = [Line2D([], [], color=rw.C_CS_TARGET, lw=LW, label=cs_label()),
               Line2D([], [], color=rw.C_DUNE_TARGET, lw=LW,
                      label="Dune-line target (net change, projected to 14 yr)"),
               Line2D([], [], color=INK, lw=LW_MODEL,
                      label=f"CASCADE ({MODEL_LABEL[key]})")] + beach_width_handles()
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False)
    sk = skill_df[skill_df["model_ends"] == folder]
    def _s(w, t):
        r = sk[(sk["window"] == w) & (sk["target"] == t)].iloc[0]
        return f"bias {r['bias_m']:+.1f} m, RMSE {r['rmse_m']:.1f} m"
    stats = " ".join(
        f"{w.replace('_', '–')}: against CoastSat {_s(w, 'coastsat')}; against the "
        f"dune line {_s(w, 'duneline')}." for w in ("1996_2010", "2010_2024"))
    metas = {o.window: o.meta for o in observations}
    dates = "; ".join(f"{m['start_date']} to {m['end_date']} ({m['interval_yr']:.1f} yr) "
                      f"for {w[0]}–{w[1]}" for w, m in metas.items())
    caption(fig, (
        "The two candidate targets and the CASCADE hindcast as net change in "
        "shoreline position over each 14-yr model window by GIS domain (1 at Cape "
        "Point, 90 at Pea Island), seaward positive, domain means. Blue: the "
        f"CoastSat target, {cs_clause()}, multiplied by 14 yr. Red: the dune-line target, the measured net "
        f"change between the digitized lines ({dates}; the 2023 date assumed) "
        "divided by its interval and multiplied by 14 yr. Black: the model's own "
        f"net change over the window, the edgeBE run with its {MODEL_LABEL[key]}, "
        "full management, groin off. The space between the two targets is the "
        "beach-width change they imply: solid grey where the beach widened, hatched "
        "where it narrowed. Interior GIS 2–89, model minus target: " + stats
        + " The LOESS-smoothed scores, as the runs are graded, are in "
        f"tables/skill.csv. One y axis, ±{half:g} m."))
    out = save(fig, OUT_DIR / folder / "target_comparison_1996_2010_2024")
    plt.close(fig)
    return out


PAIR_KEYS = (("coastsat", "coastsat_target_m", rw.C_CS_TARGET,
              None,
              "CASCADE, ends solved on CoastSat"),
             (rw.MAIN_DUNE, "dune_target_m", rw.C_DUNE_TARGET,
              "Dune-line target (net change, projected to 14 yr)",
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


def paired_figure(observations, frames, half, skill_df, ends):
    """Each target with the run solved on it, ONE FIGURE PER WINDOW, one panel
    per target (Hannah, 2026-09-19, style B of three rendered candidates):
    the target as the house fill (blue seaward, red landward), its run as the
    black line, the misfit the gap between them. Both windows on one y axis."""
    for (key, w), (_, _, n) in ends.items():
        if n != 2:
            raise SystemExit(f"{key} {w}: {n} nonzero source/sink domains, expected "
                             "the two ends only; the caption would be wrong")
    sk = skill_df.set_index(["window", "model_ends", "target"])
    names = {"coastsat": ("CoastSat target", "coastsat", "ends_solved_on_coastsat"),
             rw.MAIN_DUNE: ("Dune-line target", "duneline", "ends_solved_on_duneline")}
    out = []
    for o in observations:
        w = "{}_{}".format(*o.window)
        df = frames[o.window]
        x = df["domain_number"].to_numpy(float)
        fig, axes = plt.subplots(2, 1, sharex=True, sharey=True, constrained_layout=True,
                                 figsize=figsize("double", height=5.6))
        for i, (ax, (key, col, _, _, _)) in enumerate(zip(axes, PAIR_KEYS)):
            obs.draw_panel(ax, df.assign(mean_lrr=df[col], std_lrr=0.0), half,
                           label=(i == 0), std=False)
            obs.draw_shoals(ax, label=(i == 0))
            ax.plot(x, df[f"model_{MODEL_SETS[key]}_m"], color=INK, lw=LW_MODEL, zorder=12)
            ax.yaxis.set_major_locator(MultipleLocator(20.0 if half > 60 else 10.0))
            _title(ax, i, f"{names[key][0]} and its run")
        fills = obs.fills_in(*o.window)
        if fills:
            obs.draw_fills(axes[0], fills, half)
        axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)
        fig.supylabel(Y_LABEL, fontsize=9)
        fig.legend(handles=[(Line2D([], [], color=obs.C_ACCRETE, lw=1.0),
                             Line2D([], [], color=obs.C_ERODE, lw=1.0)),
                            Line2D([], [], color=INK, lw=LW_MODEL)],
                   labels=["Target, net change over 14 yr (seaward / landward)",
                           "CASCADE, ends solved on that target"],
                   loc="outside lower center", ncol=2, frameon=False,
                   handler_map={tuple: HandlerTuple(ndivide=None, pad=0.3)})
        (c1, c90, _), (d1, d90, _) = ends[("coastsat", o.window)], ends[(rw.MAIN_DUNE, o.window)]
        fig.suptitle("{}–{}: source/sink correction at the end domains (GIS 1 and 90) only; "
                     "none on GIS 2–89".format(*o.window) + chr(10)
                     + f"end terms GIS 1 / 90 (m/yr): (a) {c1:+.1f} / {c90:+.1f}   "
                       f"(b) {d1:+.1f} / {d90:+.1f}"
                     + ("   ·   CoastSat target = 1996–2024 LRR × 14 yr"
                        if CS_MODE == "full" else ""), fontsize=9)
        m = o.meta
        caption(fig, (
            f"{o.window[0]}–{o.window[1]}: each candidate target with the CASCADE run "
            "calibrated to it, as net change in shoreline position over the 14-yr model "
            "window by GIS domain (1 at Cape Point, 90 at Pea Island), seaward positive, "
            "domain means. (a) The CoastSat target, " + cs_clause() + " "
            "x 14 yr, as the fill (blue seaward, red landward), and the edgeBE run whose "
            "two end domains were solved against it (black). (b) The dune-line target, "
            f"the measured net change between the digitized lines ({m['start_date']} to "
            f"{m['end_date']}, {m['interval_yr']:.1f} yr"
            + (", the end date assumed" if m['end_date_assumed'] else "")
            + ") divided by its interval and x 14 yr, as the fill, and the edgeBE run whose "
            "end domains were solved against the dune line (black). The model lines are "
            "each run's own net change, unchanged; the gap between line and fill is the "
            "misfit. THE ONLY SOURCE/SINK CORRECTION IN EITHER RUN IS THE BOUNDARY TERM AT "
            "GIS 1 AND GIS 90, whose values are in the figure title; every domain from GIS "
            "2 to 89 carries none, so the interior is the model's own response. Full "
            "management, groin off. Interior GIS 2–89, model minus its own target: (a) "
            "{:+.1f} m bias, {:.1f} m RMSE; (b) {:+.1f} m bias, {:.1f} m RMSE. The y axis "
            "(±{:g} m) is shared with the other window's figure. Scores against the other "
            "target and the LOESS-smoothed targets are in tables/skill.csv.".format(
                sk.loc[(w, "ends_solved_on_coastsat", "coastsat"), "bias_m"],
                sk.loc[(w, "ends_solved_on_coastsat", "coastsat"), "rmse_m"],
                sk.loc[(w, "ends_solved_on_duneline", "duneline"), "bias_m"],
                sk.loc[(w, "ends_solved_on_duneline", "duneline"), "rmse_m"], half)))
        out += save(fig, OUT_DIR / "paired" / f"target_and_own_run_{w}")
        plt.close(fig)
    return out


CS_SOURCE = {}   # window -> (domain frame, LOESS target frame) when full-period


def main() -> int:
    global CS_MODE, OUT_DIR
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--coastsat-target", choices=sorted(CS_MODES), default="full")
    CS_MODE = ap.parse_args().coastsat_target
    OUT_DIR = ROOT_DIR / CS_MODES[CS_MODE]
    apply_style()
    observations = [rw.Observation(w) for w in WINDOWS]
    models = {}
    for key in MODEL_SETS:
        if key == "coastsat" and CS_MODE == "full":
            runs = full_period_runs()
            loaded = [rw.load_model(w, runs[w], key) for w in rw.WINDOWS]
            mdfs, rows = [m for m, _ in loaded], [r for _, r in loaded]
        else:
            mdfs, rows = rw.load_models(key)
        models[key] = ({w: m for w, m in zip(rw.WINDOWS, mdfs)},
                       [r for r in rows if tuple(int(x) for x in r["window"].split("_")) in WINDOWS])
    if CS_MODE == "full":
        full = (obs.load_window(*FULL_WINDOW), rw.load_coastsat_target(FULL_WINDOW))
        CS_SOURCE.update({o.window: full for o in observations})
    frames, skill = {}, []
    for o in observations:
        df = window_values(o, {k: m[0][o.window] for k, m in models.items()})
        frames[o.window] = df
        skill += skill_rows(o.window, df)
    skill_df = pd.DataFrame(skill)

    tables = OUT_DIR / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    for w, df in frames.items():
        df.round(3).to_csv(tables / "domain_values_{}_{}.csv".format(*w), index=False)
    skill_df.to_csv(tables / "skill.csv", index=False)
    pd.DataFrame([dict(r, model_ends=MODEL_SETS[k]) for k, (_, rows) in models.items()
                  for r in rows]).to_csv(OUT_DIR / "runs_used.csv", index=False)

    cols = ["coastsat_target_m", "dune_target_m"] + [f"model_{f}_m" for f in MODEL_SETS.values()]
    ext = max(float(np.nanmax(np.abs(df[cols].to_numpy(float)))) for df in frames.values())
    half = float(math.ceil((ext + 5) / 10.0) * 10.0)
    written = []
    for key, folder in MODEL_SETS.items():
        written += figure(observations, frames, key, folder, half, skill_df)
    written += paired_figure(observations, frames, half, skill_df,
                             end_values({k: rows for k, (_, rows) in models.items()}))

    print(skill_df[skill_df["target"].isin(["coastsat", "duneline"])].to_string(index=False))
    print(f"\ny axis +/-{half:g} m")
    for p in written:
        print("wrote   ", Path(p).relative_to(_REPO))
    return 0


if __name__ == "__main__":
    sys.exit(main())
