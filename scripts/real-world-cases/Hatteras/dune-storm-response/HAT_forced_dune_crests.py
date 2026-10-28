"""
Does the missed erosion north of Buxton come from dunes that never came down after Dennis?

    python scripts/real-world-cases/Hatteras/dune-storm-response/HAT_forced_dune_crests.py run
    python scripts/real-world-cases/Hatteras/dune-storm-response/HAT_forced_dune_crests.py score
    python scripts/real-world-cases/Hatteras/dune-storm-response/HAT_forced_dune_crests.py plot
    python scripts/real-world-cases/Hatteras/dune-storm-response/HAT_forced_dune_crests.py crests

Barrier3D lowers a dune cell only when Rhigh overtops it, so Dennis (Rhigh 3.65 m NAVD88) leaves the
5-7 m modelled dunes at GIS 8-12 standing, while the 1999 post-Dennis lidar has them at 1.7-4.7 m.
This overwrites the modelled crests with lidar survey crests (survey-crests/survey_dune_crests.csv) in
the survey's model state, in-process only, then lets the unchanged hindcast carry on to 2009: the 1998
survey (pre-Dennis) into 1 Jan 1999, the 1999 survey (post-Dennis, before Floyd) into 1 Jan 2000.
Each cell's two dune rows are scaled together so the higher one equals the survey crest; cells
without a good survey crest keep the model's. Managed edgeBE 1996-2009 runs with the groin, as the
headline matrix run:

    control       nothing forced (the matrix run, re-run here so the code is the same)
    buxton        GIS 8-12, 1999 survey
    island        every surveyed cell, 1999 survey
    buxton_pre    GIS 8-12, 1998 survey: the dunes Dennis met, so Dennis overwashes where they were low
    buxton_both   GIS 8-12, 1998 then 1999 survey
    buxton_8_15   GIS 8-15, 1999 survey: the reach CoastSat erodes, past GIS 12 where the 8-12 runs stop
    events        every surveyed cell, at the key events: 1998 (pre-Dennis), 2000 (after Dennis and
                  Floyd, the state the 1999 survey was forced into), 2004 (after Isabel)
    buxton_events GIS 8-15, the same three surveys
    all_surveys   every surveyed cell, every survey 1997-2005: the dunes as observed

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-26
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
# Hatteras/ holds the scripts; data/ and output/ live in the workspace on D: (hatteras_paths)
_HAT_DIR = next(_p for _p in _HERE.parents if (_p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in sys.path:
    sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE as PROJECT_ROOT  # noqa: E402
from cascade_pipeline.saved_runs import load_saved_cascade  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
TAG = "dunes/2026-10-26-forced-dune-crests"
EXP_DIR = PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
HINDCAST = _HAT_DIR / "hatteras_ms" / "HAT_hindcast.py"
SURVEY_CRESTS = (PROJECT_ROOT / "data" / "hatteras_init" / "1-barrier3d-domains" / "1996-start"
                 / "survey-crests" / "survey_dune_crests.csv")

START, END, PRESET, SCENARIO = 1996, 2009, "edgeBE", "full_management"
# Survey -> updates run before it is forced: the survey's model_state
#   1997 (Sept-Oct 1997)                          after 1997, the state of 1 Jan 1998
#   1998 (Sept 1998, after Bonnie, pre-Dennis)    after 1997-1998, the state of 1 Jan 1999
#   1999 (Sept 1999, post-Dennis, before Floyd)   after 1997-1999, the state of 1 Jan 2000
#   2000 (Feb-Mar 2000, after Dennis and Floyd)   after 1997-1999, the state of 1 Jan 2000
#   2001 (Jan-Mar 2001)                           after 1997-2000, the state of 1 Jan 2001
#   2004 (Apr-Sep 2004, after Isabel)             after 1997-2003, the state of 1 Jan 2004
#   2005 (Aug-Nov 2005, after Ophelia)            after 1997-2005, the state of 1 Jan 2006
FORCE_AFTER_UPDATES = {1997: 1, 1998: 2, 1999: 3, 2000: 3, 2001: 4, 2004: 7, 2005: 9}
MHW_NAVD88_M = 0.36
BUXTON = range(8, 13)
# GIS forced and the surveys forced into them
VARIANTS = {"control": None,
            "buxton": (BUXTON, (1999,)),
            "island": (range(1, 91), (1999,)),
            "buxton_pre": (BUXTON, (1998,)),
            "buxton_both": (BUXTON, (1998, 1999)),
            "buxton_8_15": (range(8, 16), (1999,)),
            # Key events, island-wide: the pre-Dennis dunes, then after Dennis and Floyd, then after Isabel
            "events": (range(1, 91), (1998, 2000, 2004)),
            "buxton_events": (range(8, 16), (1998, 2000, 2004)),
            # Every survey before the end: the dunes as observed, the most the dunes can explain
            "all_surveys": (range(1, 91), (1997, 1998, 2000, 2001, 2004, 2005))}
REPORT_GIS = range(5, 21)
WORKERS = 3
# -----------------------------------------------------------------------------


# Survey crests (m NAVD88) by GIS and cell, the ok cells only
def survey_crests(survey):
    s = pd.read_csv(SURVEY_CRESTS)
    s = s[(s.survey == survey) & (s.status == "ok")]
    return {g: dict(zip(grp.cell, grp.crest_m_navd88)) for g, grp in s.groupby("gis")}


# Overwrites the latest dune state of each listed domain with the survey crests
def force_crests(cascade, gis_list, dom, crests):
    dom = dom.for_total(len(cascade.barrier3d))
    n = 0
    for g in gis_list:
        if g not in crests:
            continue
        b = cascade.barrier3d[dom.gis_to_pad(g)]
        t = b._time_index - 1
        for cell, crest in crests[g].items():
            target = max((crest - MHW_NAVD88_M) / 10.0 - b._BermEl, b._DuneRestart)   # dam above the berm
            rows = b._DuneDomain[t, cell, :]
            top = rows.max()
            b._DuneDomain[t, cell, :] = (np.maximum(rows * target / top, b._DuneRestart)
                                         if top > 0 else target)
            n += 1
    return n


# Child process: force each survey's crests in its year, then run the unchanged hindcast
def _launch(variant):
    import runpy
    import cascade_pipeline.hindcast as H
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS
    if VARIANTS[variant] is not None:
        gis_list, surveys = VARIANTS[variant]
        due = {FORCE_AFTER_UPDATES[s]: (s, survey_crests(s)) for s in surveys}
        orig_build = H.build_cascade

        def build_cascade(*a, **k):
            c = orig_build(*a, **k)
            cls = type(c)
            orig_update = cls.update
            calls = [0]

            # Patched on the class, so the saved model still pickles
            def update(self):
                orig_update(self)
                calls[0] += 1
                if calls[0] in due:
                    survey, crests = due[calls[0]]
                    n = force_crests(self, gis_list, HATTERAS_DOMAINS, crests)
                    print(f"\nEXPERIMENT {n} dune cells forced to the {survey} survey crests "
                          f"(GIS {min(gis_list)}-{max(gis_list)}) after update {calls[0]}", flush=True)

            cls.update = update
            return c

        H.build_cascade = build_cascade
    sys.path.insert(0, str(HINDCAST.parent))
    sys.argv = [str(HINDCAST)]
    runpy.run_path(str(HINDCAST), run_name="__main__")


def run_dir(variant):
    base = EXP_DIR / "runs" / variant / f"{START}_{END}" / PRESET
    hits = sorted(d for d in base.glob("HAT_*") if any(d.glob("*_run_metadata.json"))) if base.exists() else []
    return hits[-1] if hits else None


# One run in its own process, logged, with its launch record
def launch(variant):
    logs = EXP_DIR / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env.update({"HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(START),
                "HAT_SOURCE_SINK_PRESET": PRESET, "HAT_SCENARIO": SCENARIO, "HAT_GROIN_ENABLED": "1",
                "HAT_RUN_KIND": "experiment", "HAT_RUN_TAG": f"{TAG}/runs/{variant}",
                "HAT_SAVE_MODEL_STATE": "1", "HAT_MAKE_GIFS": "0", "HAT_SHOW_FIGURES": "0",
                "HAT_OVERWRITE": "1", "PYTHONIOENCODING": "utf-8"})
    t0 = time.time()
    log = logs / f"{variant}.log"
    with open(log, "w", encoding="utf-8") as fh:
        p = subprocess.run([sys.executable, str(_HERE), "_launch", variant],
                           env=env, cwd=HINDCAST.parent, stdout=fh, stderr=subprocess.STDOUT)
    rec = dict(variant=variant, surveys=VARIANTS[variant] and list(VARIANTS[variant][1]), returncode=p.returncode,
               minutes=round((time.time() - t0) / 60, 1))
    with open(logs / "launches.jsonl", "a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec) + "\n")
    print(f"{variant:8s} exit {p.returncode}  {rec['minutes']} min", flush=True)


# Scores

# Per GIS: LRR, 2000-2009 endpoint rate, overwash and median crest, against the CoastSat target
def score():
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    frames = {}
    for v in VARIANTS:
        rd = run_dir(v)
        if rd is None:
            print(f"  missing {v}")
            continue
        c = load_saved_cascade(rd / f"{rd.name}.npz")
        dom = DOM.for_total(len(c.barrier3d))
        rate = pd.read_csv(rd / "tables" / "shoreline_change_rate.csv").set_index("gis_domain")
        sh = -np.load(next(rd.glob("*_shoreline_matrix.npy")))
        ow = np.load(next(rd.glob("*_overwash.npy")))
        meta = json.loads((rd / f"{rd.name}_run_metadata.json").read_text(encoding="utf-8"))
        rows = []
        for g in range(1, 91):
            p = dom.gis_to_pad(g)
            b = c.barrier3d[p]
            crest = (np.asarray(b.DuneDomain, float).max(axis=2) + b._BermEl) * 10 + MHW_NAVD88_M
            rows.append(dict(gis=g, lrr=rate.lrr_m_yr[g],
                             rate_2000_2009=(sh[-1, p] - sh[3, p]) / (len(sh) - 1 - 3),
                             overwash_2000_2009=ow[4:, p].sum(), crest_2000=np.median(crest[3]),
                             crest_2009=np.median(crest[-1])))
        frames[v] = pd.DataFrame(rows).set_index("gis")
        print(f"{v:12s} interior bias {float(meta['skill']['mean_bias_interior_m_yr']):+.3f}  "
              f"rmse {float(meta['skill']['rmse_interior_m_yr']):.3f} m/yr")
    if not frames:
        return
    out =pd.concat(frames, axis=1).loc[list(REPORT_GIS)]
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 40)
    print(out.round(2))
    out.to_csv(EXP_DIR / "forced_dune_crests_by_gis.csv")
    print(f"\nwrote {EXP_DIR / 'forced_dune_crests_by_gis.csv'}")


# Figures

# Colour and line width per run, all solid: the two forced runs one green, light to dark, as the same
# 1999 dunes over a wider and a narrower reach; drawn in this order, so later runs sit on top where they coincide
RUN_STYLES = {
    "control": ("No forcing (control)", "#8C8C8C", 1.6),
    "island": ("1999 lidar dunes, all GIS", "#74C476", 1.8),
    "buxton_8_15": ("1999 lidar dunes, GIS 8-15", "#006D2C", 2.2),
}
# (file tag, GIS range or None for the whole island, y half-range in m)
VIEWS = (("island", None, 100.0), ("buxton", (1, 25), 60.0))
FORCED_GIS = (8, 15)
# Map features with nothing to point at in the Buxton view
ZOOM_HIDDEN_FEATURES = {"Shoal influence zone", "Village center", "Pier"}


# Every run's end-minus-start position change against the CoastSat net change, in the layout of the run's
# shoreline_position_change figure (net_change_comparison.py); one figure per view in VIEWS
def plot():
    import dataclasses
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from cascade_pipeline.coastsat_lowess import CoastSatDataset, build_coastsat_series
    from cascade_pipeline.plotting.comparison_style import restyle
    from cascade_pipeline.plotting.net_change_comparison import LOWESS_CONFIG, smooth_model
    from cascade_pipeline.plotting.rate_comparison import (
        DEFAULT_RATE_COMPARISON, plot_annotated_rate_comparison)
    from cascade_pipeline.run_info import RunInfo
    from cascade_pipeline.run_layout import resolve
    from cascade_pipeline.shoreline import compute_change_rate
    from site_layer.hat_figure_style import mark_offaxis, offaxis_clause, record_caption
    from site_layer.hat_observed_rates import NET_CHANGE_WINDOWS, WINDOW_ROLE, net_change_dir
    from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS, HATTERAS_DOMAINS

    # End minus start (m, seaward positive), smoothed as the run figure smooths it
    runs = {}
    for v in RUN_STYLES:
        rd = run_dir(v)
        if rd is None:
            print(f"  missing {v}")
            continue
        shoreline_m = np.load(resolve(rd, "matrix", rd.name))
        runs[v] = (rd, smooth_model(compute_change_rate(shoreline_m, span_years=1, flip_sign=True)))
    if "control" not in runs:
        print("  the control run is needed to draw the figures")
        return

    tag = f"{START}_{END}"
    cs_series = build_coastsat_series(
        [CoastSatDataset(label=f"CoastSat net change ({START}–{END})", period_start=START,
                         csv_path=str(net_change_dir(START, END) / f"transect_net_change_{tag}.csv"),
                         rate_col="net_change_m")],
        active_period_start=START, lowess_config=LOWESS_CONFIG, domains=HATTERAS_DOMAINS)

    control_dir, control = runs["control"]
    label_c, color_c, _ = RUN_STYLES["control"]
    gis = np.arange(HATTERAS_DOMAINS.first_gis_id, HATTERAS_DOMAINS.last_gis_id + 1)
    geometry = HATTERAS_DOMAINS.for_total(len(control))     # the runs' own buffers, not the config's
    real = slice(geometry.start_real_index, geometry.end_real_index)
    role = WINDOW_ROLE.get((START, END), "").replace(" (held out)", "")

    for view, gis_range, y_half in VIEWS:
        where = f"GIS {gis_range[0]}–{gis_range[1]}" if gis_range else "whole island"
        subtitle = " · ".join(p for p in (role, where, "full management + groin, edgeBE",
                                          "dune crests forced to the 1999 post-Dennis lidar") if p)
        config = dataclasses.replace(
            DEFAULT_RATE_COMPARISON, quantity="position", publication_text=True,
            show_features=True, plot_domain_means=False, raw_lrr_southern_only=True,
            smoothing_marker=False, observed_label="CoastSat (7-domain LOWESS)", model_label=label_c,
            title=f"Shoreline position change with the dunes forced to the lidar, {START}–{END}",
            ylim=(-y_half, y_half), subtitle=subtitle)
        run = RunInfo(run_name=control_dir.name, run_dir=str(control_dir), start_year=START,
                      end_year=END, Hs=2.0, flip_sign_model=True)
        # The control is drawn by the run figure's own call, in the control's colour
        fig, ax = plot_annotated_rate_comparison(
            control, cs_series, run, estimator=None, domains=HATTERAS_DOMAINS,
            annotations=dataclasses.replace(HATTERAS_ANNOTATIONS, model_color=color_c),
            lowess_config=LOWESS_CONFIG, config=config)

        in_view = np.ones_like(gis, bool) if gis_range is None else (
            (gis >= gis_range[0]) & (gis <= gis_range[1]))
        named = [(label_c, mark_offaxis(ax, gis[in_view], control[real][in_view], y_half, color=color_c))]
        handles = []
        for v, (rd, change) in runs.items():
            label, color, lw = RUN_STYLES[v]
            handles.append(Line2D([0], [0], color=color, lw=lw + 0.4, label=label))
            if v == "control":
                # Drawn by the run figure's call; given its width here, under the forced runs
                for line in ax.get_lines():
                    if line.get_label() == label_c:
                        line.set_linewidth(lw)
                continue
            ax.plot(gis, change[real], color=color, lw=lw, zorder=6)
            named.append((label, mark_offaxis(ax, gis[in_view], change[real][in_view], y_half, color=color)))
        ax.axvspan(FORCED_GIS[0] - 0.5, FORCED_GIS[1] + 0.5, color=RUN_STYLES["buxton_8_15"][1],
                   alpha=0.07, lw=0, zorder=0)

        if gis_range is not None:
            # Zoom: hide the labels of features outside the view; centre the transects on their
            # domain, as groin_effect_buxton.py does
            ax.set_xlim(gis_range[0] - 0.5, gis_range[1] + 0.5)
            ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
            ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
            for txt in ax.texts:
                x_is_data = txt.get_transform() not in (ax.transAxes, fig.transFigure)
                if x_is_data and not gis_range[0] - 0.5 <= txt.get_position()[0] <= gis_range[1] + 0.5:
                    txt.set_visible(False)
            for coll in ax.collections:
                if "transect" in str(coll.get_label()):
                    coll.set_offsets(coll.get_offsets() - [0.5, 0.0])

        # The figure's legend with every run in place of the one model line, then the house finish
        old = fig.legends.pop()
        hidden = {label_c} | (set() if gis_range is None else ZOOM_HIDDEN_FEATURES)
        entries = [(h, t.get_text()) for h, t in zip(old.legend_handles, old.get_texts())
                   if t.get_text() not in hidden]
        handles.append(Patch(fc=RUN_STYLES["buxton_8_15"][1], alpha=0.15,
                             label=f"Dunes forced (GIS {FORCED_GIS[0]}-{FORCED_GIS[1]} run)"))
        fig.legend(handles + [h for h, _ in entries],
                   [h.get_label() for h in handles] + [t for _, t in entries],
                   loc="outside lower center")
        restyle(fig, ax, subtitle)

        out = EXP_DIR / f"forced_dune_crests_position_{view}.png"
        fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
        plt.close("all")
        (s0, s1), (e0, e1) = NET_CHANGE_WINDOWS[(START, END)]
        record_caption(out, (
            f"Shoreline position change, {START}–{END}, {where}; positive is seaward. Modelled "
            f"(CASCADE): the shoreline at the end of each run minus the start, smoothed with a "
            f"7-domain LOWESS north of domain 10 as in the run figure, for the full-management edgeBE "
            f"run with the groin and nothing forced (grey) and with the dune crests overwritten by the "
            f"1999 post-Dennis lidar survey crests in the state of 1 Jan 2000, at GIS "
            f"{FORCED_GIS[0]}-{FORCED_GIS[1]} (dark green) or at every surveyed domain (light green). Green "
            f"shading: GIS {FORCED_GIS[0]}-{FORCED_GIS[1]}. Where runs coincide "
            f"the last drawn sits on top. Observed (CoastSat, blue): per transect, the mean position "
            f"over {e0} to {e1} minus the mean over {s0} to {s1}, smoothed with a 7-domain LOWESS "
            f"north of domain 10; dots south of it. The y axis is fixed at ±{y_half:g} m."
            + offaxis_clause([(k, n) for k, n in named if n], y_half)
            + f" Experiment {TAG}."))
        print(f"wrote {out}")


# When the modelled dunes leave the lidar: the control's median crest each year against every survey's,
# one panel per domain of the forced reach

CREST_GIS = range(8, 16)
STORMS = {"Bonnie": 1998.65, "Dennis": 1999.67, "Isabel": 2003.71}   # landfall, decimal year


def crests():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import FixedLocator, FormatStrFormatter, MultipleLocator
    from site_layer.hat_figure_style import apply_style, record_caption
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS
    apply_style()

    rd = run_dir("control")
    c = load_saved_cascade(rd / f"{rd.name}.npz")
    dom = HATTERAS_DOMAINS.for_total(len(c.barrier3d))
    s = pd.read_csv(SURVEY_CRESTS)
    s = s[s.status == "ok"].copy()
    s["row"] = pd.to_datetime(s.model_state).dt.year - (START + 1)    # row 0 is 1 Jan of the first model year
    # Mid-point of the months flown, as a decimal year
    span = s.flown.str.split(" to ")
    first, last = (pd.to_datetime(span.str[0]), pd.to_datetime(span.str[-1]))
    s["when"] = (first.dt.year + (first.dt.month + last.dt.month) / 2 / 12
                 + (last.dt.year - first.dt.year) * last.dt.month / 24)

    model_c, lidar_c = RUN_STYLES["control"][1], "#2166AC"
    fig, axs = plt.subplots(2, 4, figsize=(12, 6.2), sharex=True, sharey=True)
    rows = []
    for ax, g in zip(axs.flat, CREST_GIS):
        b = c.barrier3d[dom.gis_to_pad(g)]
        crest = (np.asarray(b.DuneDomain, float).max(axis=2) + b._BermEl) * 10 + MHW_NAVD88_M
        years = START + 1 + np.arange(crest.shape[0])
        lo, mid, hi = np.percentile(crest, [25, 50, 75], axis=1)
        ax.fill_between(years, lo, hi, color=model_c, alpha=0.18, lw=0)
        ax.plot(years, mid, color=model_c, lw=1.8)
        for sv, grp in s[s.gis == g].groupby("survey"):
            q = grp.crest_m_navd88.quantile([0.25, 0.5, 0.75]).values
            ax.errorbar(grp.when.iloc[0], q[1], yerr=[[q[1] - q[0]], [q[2] - q[1]]], fmt="o", ms=4,
                        color=lidar_c, elinewidth=1, capsize=0, zorder=5)
            r = int(grp.row.iloc[0])
            rows.append(dict(gis=g, survey=sv, lidar=q[1],
                             model=float(np.median(crest[r][grp.cell.values]))))
        for name, t in STORMS.items():
            ax.axvline(t, color="0.55", lw=0.8, ls=":")
            if g == CREST_GIS[0]:
                ax.text(t, 0.03, name, rotation=90, ha="right", va="bottom", fontsize=7, color="0.4",
                        transform=ax.get_xaxis_transform())
        ax.set_title(f"GIS {g}", loc="left", fontsize=10, fontweight="bold")
    for ax in axs[:, 0]:
        ax.set_ylabel("Dune crest (m NAVD88)")
    for ax in axs[-1]:
        ax.set_xlabel("Year")
    # Whole years: a labelled tick on every even year, a minor tick each year
    axs[0, 0].set_xlim(START - 0.5, END + 1.5)
    axs[0, 0].xaxis.set_major_locator(FixedLocator(range(START, END + 2, 2)))
    axs[0, 0].xaxis.set_minor_locator(MultipleLocator(1))
    axs[0, 0].xaxis.set_major_formatter(FormatStrFormatter("%d"))
    fig.legend(handles=[Line2D([], [], color=model_c, lw=2, label="Model, control (median; band: 25th-75th pct)"),
                        Line2D([], [], color=lidar_c, marker="o", ls="none",
                               label="Lidar survey (median; bar: 25th-75th pct)")],
               loc="lower center", ncol=2, frameon=False)
    fig.suptitle(f"Modelled and surveyed dune crests north of Buxton, {START}–{END}",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    out = EXP_DIR / "dune_crests_model_vs_lidar_gis8_15.png"
    fig.savefig(out, dpi=250, facecolor="white")
    plt.close(fig)
    record_caption(out, (
        f"Dune crest elevation (m NAVD88) per GIS domain, {CREST_GIS[0]}-{CREST_GIS[-1]}: the control "
        f"run's annual state (1 Jan; line the median over the domain's 50 cells, band the 25th-75th "
        f"percentile) and each lidar survey's crests on the same cells (survey-crests/survey_dune_crests.csv; "
        f"dot the median, bar the 25th-75th percentile), at the middle of the months flown. Dotted: "
        f"Bonnie, Dennis and Isabel. 2001 is unclassified lidar, so buildings or vegetation can raise "
        f"its crests. Run {rd.name}, experiment {TAG}."))
    d = pd.DataFrame(rows)
    d["model_minus_lidar"] = d.model - d.lidar
    d.round(3).to_csv(EXP_DIR / "dune_crests_model_vs_lidar_gis8_15.csv", index=False)
    print(d.pivot(index="gis", columns="survey", values="model_minus_lidar").round(1))
    print(f"wrote {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("action", choices=["run", "score", "plot", "crests", "_launch"])
    ap.add_argument("variant", nargs="*")
    args = ap.parse_args()
    if args.action == "_launch":
        _launch(args.variant[0])
    elif args.action == "run":
        variants = args.variant or list(VARIANTS)
        with ThreadPoolExecutor(WORKERS) as ex:
            list(ex.map(launch, variants))
    elif args.action == "plot":
        plot()
    elif args.action == "crests":
        crests()
    else:
        score()


if __name__ == "__main__":
    main()
