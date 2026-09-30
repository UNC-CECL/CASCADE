"""Dune line vs shoreline as the island offset (orientation), 1996-2010 (2026-09-25).

Asked by Hannah on 2026-09-25: how much does setting the island's planform
from the dune line rather than the CoastSat shoreline change the output?

    offsets   duneline (1996/duneline/v1) and shoreline (1996/shoreline/v1),
              both metres, the Hermite wrap-around written in the file
    waves     Hs 1.0 m, Tp 8 s, asymmetry 0.8 -- the best managed 1996-2010
              setting found (tuned with the DUNE-LINE offset) -- and, so the
              shoreline offset gets a fair chance, high-angle fraction
              0.3 0.4 0.45 0.5 0.55 (the lever that mattered most)
    scope     natural and full management, 1996-2010
    score     as wave-climate/2026-09-25-wave-grid-smoothed-score: share of the alongshore
              variation explained by the model SMOOTHED like the CoastSat
              target, interior GIS 2-89; raw score, bias and r beside it
Both offsets are run fresh here (20 runs) so every run in the comparison is
on the same code and the same Barrier3D (the route_overwash fix).

WHERE: output/raw_runs/experiments/island-offset/2026-09-25-metres-offset-duneline-vs-shoreline-waves-hs1-tp8/
    README.md, tables/all_runs.csv, figures/, logs/<source>_<scenario>/<settings>.log
    runs/<source>_<scenario>/1996_2010/zeroBE/<run_name>/   (on disk only)

    python scripts/hatteras_ms/experiments/HAT_offset_source_comparison.py run
    python scripts/hatteras_ms/experiments/HAT_offset_source_comparison.py score
    python scripts/hatteras_ms/experiments/HAT_offset_source_comparison.py plot
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from itertools import product
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_wave_grid_smoothed_score as grid  # noqa: E402

common, step2 = grid.common, grid.step2
TAG = "island-offset/2026-09-25-metres-offset-duneline-vs-shoreline-waves-hs1-tp8"
STUDY_DIR = grid.RAW_RUNS / "experiments" / TAG
TABLES_DIR, LOGS_DIR, FIG = STUDY_DIR / "tables", STUDY_DIR / "logs", STUDY_DIR / "figures"
PERIOD = 1996
SOURCES = ("duneline", "shoreline")
SCENARIOS = grid.SCENARIOS
BASE = {"hs": 1.0, "wave_period_s": 8.0, "wave_asymmetry": 0.8}
HIGH_ANGLE = (0.3, 0.4, 0.45, 0.5, 0.55)
HEADLINE = 0.45
# Figure labels, overridden by a study that reuses this driver
# (HAT_offset_source_comparison_div10.py).
WAVE_NOTE = ("the 09-25 setting, tuned with the dune-line offset; superseded by "
             "option A on 09-27")
OFFSET_NOTE = "in metres (duneline/v1, shoreline/v1)"
FORM_NOTE = ("Redrawn 2026-09-28 in the form of "
             "../2026-09-28-metres-offset-duneline-vs-shoreline-waves-option-a/.")
LEGEND_TITLE = "No source/sink correction at the ends"


def cells():
    return [(src, sc, {**BASE, "wave_angle_high_fraction": f})
            for src, sc, f in product(SOURCES, SCENARIOS, HIGH_ANGLE)]


def group(src, sc):
    return f"{src}_{sc}"


def log_path(src, sc, s):
    return LOGS_DIR / group(src, sc) / f"{grid.label(s)}.log"


def env(src, sc, s):
    e = grid.run_env("x", sc, PERIOD, s)
    e["HAT_ISLAND_OFFSET_SOURCE"] = src
    e["HAT_RUN_TAG"] = f"{TAG}/runs/{group(src, sc)}"
    return e


def launch(cell):
    src, sc, s = cell
    log = log_path(src, sc, s)
    if grid.finished(log):
        return
    log.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    p = subprocess.run([sys.executable, str(grid.HINDCAST)], env=env(src, sc, s),
                       cwd=str(grid.PROJECT_ROOT), capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=grid.RUN_TIMEOUT_S)
    log.write_text((p.stdout or "") + "\n--- STDERR ---\n" + (p.stderr or ""), encoding="utf-8")
    what = "done" if p.returncode == 0 else f"FAILED ({step2.stop_reason(log)})"
    print(f"{what} {group(src, sc)} {grid.label(s)} in {(time.perf_counter() - t0) / 60:.1f} min",
          flush=True)


def cmd_run(a):
    grid.check_barrier3d()
    common.keep_awake()
    todo = [c for c in cells() if not grid.finished(log_path(*c))]
    print(f"{len(cells())} cells, {len(todo)} to run, {a.jobs} at a time", flush=True)
    with ThreadPoolExecutor(max_workers=a.jobs) as pool:
        list(pool.map(launch, todo))
    return cmd_score(a)


def cmd_score(_=None):
    import pandas as pd
    from cascade_pipeline.run_registry import load_run_index, rebuild_run_index
    rebuild_run_index(grid.RAW_RUNS)
    idx = load_run_index(grid.RAW_RUNS / "run_index.csv")
    idx = idx[idx["tag"].astype(str).str.startswith(TAG + "/") & (idx["status"] == "current")]
    target = common.coastsat_target(PERIOD)
    runs = {}
    for _, r in idx.iterrows():
        d = (grid.RAW_RUNS / "experiments" / r["tag"] / f"{r['start_year']}_{r['end_year']}"
             / r["source_sink_preset"] / r["run_name"])
        md = json.loads((d / f"{r['run_name']}_run_metadata.json").read_text(encoding="utf-8"))
        runs[(r["tag"].split("/")[-1], float(md["wave climate"]["wave_angle_high_frac"]))] = (d, md)
    rows = []
    for src, sc, s in cells():
        log = log_path(src, sc, s)
        rec = {"source": src, "scenario": sc, **s}
        hit = runs.get((group(src, sc), s["wave_angle_high_fraction"]))
        if hit and grid.finished(log):
            d, md = hit
            got = md["identity"]["island_offset_version"]
            # a superseded build records the source alone
            # (hatteras_site_config.island_offset_version)
            if not (str(got) == src or str(got).startswith(src + "/")):
                raise ValueError(f"{d}: ran on offset {got!r}, filed as {src}")
            rec.update(status="scored", **grid.score_run(d, target),
                       offset_recorded=str(got),
                       run_dir=str(d.relative_to(STUDY_DIR)).replace("\\", "/"))
        else:
            rec["status"] = step2.stop_reason(log) if log.is_file() else "not run"
        rows.append(rec)
    t = pd.DataFrame(rows)
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    t.to_csv(TABLES_DIR / "all_runs.csv", index=False)
    cols = ["source", "scenario", "wave_angle_high_fraction", "status",
            "smoothed_variance_explained", "raw_variance_explained", "bias_m_yr", "smoothed_r"]
    print(t[[c for c in cols if c in t]].round(3).to_string(index=False))
    return 0


def house_figures(panels, fig_dir, legend_title, note, suffix,
                  shoreline_window="over 1995-1997 for 1996, 2009-2011 for 2010"):
    """The four island-offset figures, in ONE form for every study that asks the
    dune-line-or-shoreline question (Hannah, 2026-09-28: "ensure the figures
    among these experiments are consistent ... so it is easier to compare").

    panels   [{"label": "Natural 1996–2010", "start": 1996, "end": 2010,
               "rates": {"duneline": df, "shoreline": df}}, ...], one row each;
               df is a run's tables/shoreline_change_rate.csv indexed by
               gis_domain, or None where that arm has no run
    note     the study's own sentence for every caption (waves, offset, ends)
    suffix   the stem ending, e.g. "1996_2010" or "full_management"
    shoreline_window  the mean-shoreline windows the shoreline offset was built
             on, for the captions; the default is the v1 (calendar) builds

    Net change in metres; observations smoothed with LOESS over LOESS_DOMAINS
    (southern SKIP_SOUTHERN raw), the model unsmoothed; the model's ENABLED
    fills marked above each panel; no scores on the figures."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from site_layer.hat_figure_style import (INK, INK_MUTED, DOMAIN_AXIS_LABEL, C,
                                             apply_style, open_frame, record_caption, save,
                                             structures, support_dir, town_bands)
    from site_layer.hatteras_site_config import (HATTERAS_ANNOTATIONS,
                                                 HATTERAS_NOURISHMENT_PROJECTS)
    apply_style()
    import numpy as np
    from matplotlib.ticker import MultipleLocator
    col = {"duneline": "#1b7f6b", "shoreline": "#6a3d9a"}
    # Text sized for reading the figures side by side (Hannah, 2026-09-29:
    # "make all the text larger"); the canvas stays 16 in wide.
    rc = {"font.size": 17, "axes.titlesize": 18, "axes.labelsize": 17,
          "xtick.labelsize": 15, "ytick.labelsize": 15, "legend.fontsize": 15}
    # One y-axis for every island-offset study (09-22, 09-25, 09-28 option A),
    # set from the largest of them: the change figures, and the difference one
    CHANGE_YLIM, DIFF_YLIM = (-80, 100), (-20, 30)
    shoal_c = C["ADDED"]
    n = len(panels)

    def fills(ax, start, end):
        trans = ax.get_xaxis_transform()
        for p in sorted(HATTERAS_NOURISHMENT_PROJECTS, key=lambda p: p.year):
            if not (p.enabled and start <= p.year <= end):
                continue
            lo, hi = min(p.gis_domains), max(p.gis_domains)
            ax.plot([lo - 0.45, hi + 0.45], [1.025, 1.025], color=INK, lw=2.2,
                    solid_capstyle="butt", zorder=6, clip_on=False, transform=trans)
            ax.text((lo + hi) / 2, 1.045, f"{p.year} fill", ha="center", va="bottom",
                    fontsize=14, color=INK, zorder=6, clip_on=False, transform=trans)

    def dress(ax, j, pan, text):
        ax.axhline(0, color=INK_MUTED, lw=0.8)
        ax.set_xlim(1, 90)
        ax.grid(axis="y")
        open_frame(ax)
        town_bands(ax, fontsize=15, strip=0.08)
        for nm, (lo, hi) in HATTERAS_ANNOTATIONS.shoal_zones.items():
            ax.axvspan(lo - 0.5, hi + 0.5, color=shoal_c, alpha=0.12, lw=0, zorder=0.5)
            ax.text((lo + hi) / 2, 0.895, nm, transform=ax.get_xaxis_transform(),
                    ha="center", va="top", fontsize=15, color="#8a620e", zorder=1)
        structures(ax, label=True, label_pt=13)
        fills(ax, pan["start"], pan["end"])
        ax.set_title(f"({'abcdef'[j]})", loc="left", fontweight="bold", pad=40)
        ax.set_title(text, loc="center", pad=40)
        if j == n - 1:
            ax.set_xlabel(DOMAIN_AXIS_LABEL)

    def legend(f, handles):
        # legend_title (the end correction) goes to the caption, not the canvas:
        # it read as a heading for the legend entries (2026-09-29)
        handles = handles + [Patch(color=shoal_c, alpha=0.25, label="Shoals"),
                             Patch(color="0.90", label="Villages")]
        f.legend(handles=handles, loc="outside lower center", ncol=2, frameon=False)

    def axes_for(height):
        f, axes = plt.subplots(n, 1, figsize=(16, height * n + 1.8), sharex=True,
                               sharey=True, constrained_layout=True, squeeze=False)
        return f, axes[:, 0]

    def limits(values, fixed, step, what):
        """The FIXED y-limits, the same in every island-offset study so the
        figures compare across studies (Hannah, 2026-09-29), widened to `step`
        only if a study's data leaves them, with a warning."""
        v = np.concatenate([np.asarray(x, float) for x in values])
        v = v[np.isfinite(v)]
        lo, hi = fixed
        if v.min() < lo or v.max() > hi:
            lo = min(lo, np.floor(v.min() / step) * step)
            hi = max(hi, np.ceil(v.max() / step) * step)
            print(f"  WARNING: {what} data ({v.min():.1f} to {v.max():.1f} m) leave the "
                  f"shared axis {fixed}; widened to ({lo:g}, {hi:g}), no longer comparable")
        return lo, hi

    common = ((f" {legend_title}." if legend_title else "") + " Amber: Avon and Wimble Shoals; solid line: Buxton groin; dotted lines: "
              "Avon and Rodanthe piers; grey strip: villages; bars above a panel: the "
              "nourishments the model is given in that period, footprint and year. " + note)
    common += (" The three comparison figures share one y-axis; the difference figure "
               "has its own, symmetric about zero.")
    figs = [
        dict(src="duneline", kind="dune",
             ylab="Dune-line change (m)\n(+ seaward, − landward)",
             obs_label="Observed: dune-line change (7-domain LOESS)",
             model_label="Model: started from the dune line",
             stem="duneline_offset_vs_duneline_change",
             caption=("The model started from the DUNE-LINE island offset, against the dune "
                      "line's own change. Net change per domain, seaward positive: observed "
                      "(black) is the mean change between the digitised dune lines that bound "
                      "each period (1997-10 to 2009-05, 11.6 yr, for 1996-2010; 2009-05 to "
                      "2023-07, 14.1 yr, for 2010-2024), LOESS over 7 domains (the southern 10 "
                      "raw); the model (green) is unsmoothed, its endpoint change over the 14 "
                      "calendar years. The interval mismatch is not corrected.")),
        dict(src="shoreline", kind="total",
             ylab="Total shoreline change (m)\n(+ seaward, − landward)",
             obs_label="Observed: CoastSat LRR of the same period × 14 yr (7-domain LOESS)",
             model_label="Model: started from the shoreline (its LRR × 14 yr)",
             stem="shoreline_offset_vs_coastsat_total_change",
             caption=("The model started from the SHORELINE island offset (mean CoastSat "
                      f"shoreline {shoreline_window}), against total "
                      "shoreline change: each period's OWN CoastSat LRR, LOESS over 7 domains, "
                      "x 14 yr (black), not the 1996-2024 rate carried onto it; the model "
                      "(purple) is its own LRR x 14 yr.")),
        dict(src="shoreline", kind="projected",
             ylab="Projected shoreline change (m)\n(+ seaward, − landward)",
             obs_label="Observed: CoastSat LRR 1996–2024 × 14 yr (7-domain LOESS)",
             model_label="Model: started from the shoreline (its LRR × 14 yr)",
             stem="shoreline_offset_vs_coastsat_projected_change",
             caption=("The model started from the SHORELINE island offset, against PROJECTED "
                      "shoreline change: the long-term CoastSat LRR fitted on 1996-2024, LOESS "
                      "over 7 domains, x 14 yr (black), the same profile in every panel, "
                      "carried onto each period rather than fitted on it; the model (purple) "
                      "is the same run as in the total-change figure.")),
    ]
    # Every observed and modelled line first, so the three comparison figures
    # share ONE y-axis and read against each other (Hannah, 2026-09-29).
    series = {}
    for spec in figs:
        for j, pan in enumerate(panels):
            a, b = pan["start"], pan["end"]
            years = b - a
            rt = pan["rates"].get(spec["src"])
            obs = {"dune": lambda: smooth_loess(duneline_change(a, b)),
                   "total": lambda: coastsat_target_loess(a, f"{a}_{b}") * years,
                   "projected": lambda: coastsat_target_loess(1996, "1996_2024") * years,
                   }[spec["kind"]]()
            m = None
            if rt is not None:
                m = (rt.change_rate_m_yr if spec["kind"] == "dune" else rt.lrr_m_yr) * years
            series[(spec["stem"], j)] = (obs, m)
    # No headroom above the data: the one value past the label line is the
    # observed +96 m at GIS 1-2 (2010-2024), where no label sits; headroom for
    # it pushed the top to 160 m and flattened every line.
    ylim = limits([x.values for pair in series.values() for x in pair if x is not None],
                  CHANGE_YLIM, 20, "change")

    out = []
    with plt.rc_context(rc):
        for spec in figs:
            src = spec["src"]
            f, axes = axes_for(5.2)
            for j, pan in enumerate(panels):
                ax = axes[j]
                obs, m = series[(spec["stem"], j)]
                ax.plot(obs.index, obs.values, color=INK, lw=3.4, zorder=6)
                if m is not None:
                    ax.plot(m.index, m.values, color=col[src], lw=2.4, zorder=4)
                ax.set_ylim(*ylim)
                ax.yaxis.set_major_locator(MultipleLocator(20))
                a, b = pan["start"], pan["end"]
                how = {"total": f"  ·  observed: CoastSat LRR {a}–{b} × 14 yr",
                       "projected": "  ·  observed: CoastSat LRR 1996–2024 × 14 yr"
                       }.get(spec["kind"], "")
                dress(ax, j, pan, pan["label"] + how)
                ax.set_ylabel(spec["ylab"])
            legend(f, [Line2D([], [], color=INK, lw=3.4, label=spec["obs_label"]),
                       Line2D([], [], color=col[src], lw=2.4, label=spec["model_label"])])
            png = fig_dir / f"{spec['stem']}_{suffix}.png"
            save(f, png, dpi=300, close=True)
            record_caption(png, spec["caption"] + common)
            out.append(png)

        # where the two offsets disagree: total change, shoreline minus dune-line start
        f, axes = axes_for(4.4)
        rows, diffs = [], []
        for j, pan in enumerate(panels):
            ax, r = axes[j], pan["rates"]
            if r.get("shoreline") is not None and r.get("duneline") is not None:
                years = pan["end"] - pan["start"]
                d = (r["shoreline"].lrr_m_yr - r["duneline"].lrr_m_yr) * years
                diffs.append(d.values)
                ax.bar(d.index, d.values, width=0.85,
                       color=[col["shoreline"] if v > 0 else col["duneline"] for v in d.values])
                rows.append(dict(panel=pan["label"], mean_abs_diff_m=float(d.abs().mean()),
                                 max_abs_diff_m=float(d.abs().max()),
                                 at_gis=int(d.abs().idxmax())))
            dress(ax, j, pan, pan["label"])
            ax.set_ylabel("Shoreline start minus\ndune-line start (m)")
        if diffs:
            # its own quantity, so its own axis, with room above the bars for
            # the village and shoal labels
            lo, hi = limits(diffs, DIFF_YLIM, 5, "difference")
            for ax in axes:
                ax.set_ylim(lo, hi)
                ax.yaxis.set_major_locator(MultipleLocator(5))
        legend(f, [Patch(color=col["shoreline"], label="Shoreline start more accretional"),
                   Patch(color=col["duneline"], label="Dune-line start more accretional")])
        png = fig_dir / f"total_change_difference_shoreline_minus_duneline_{suffix}.png"
        save(f, png, dpi=300, close=True)
        pd.DataFrame(rows).to_csv(support_dir(png.parent) / f"{png.stem}.csv", index=False)
        record_caption(png, (
            "Modelled total shoreline change (each run's own LRR x 14 yr) with the shoreline "
            "start minus that with the dune-line start, per domain. Purple: the shoreline start "
            "gives the more accretional (less erosional) change; green: the dune-line start "
            "does. Both runs give the same model quantity, so no observation enters." + common))
        out.append(png)
    return out


# THE FIGURES ARE FULL MANAGEMENT x BOTH PERIODS (Hannah, 2026-09-28: "showing
# only full management and both periods per figure", as the option A study
# draws them). 1996's run is the study's own; 2010's comes from `run-2010`.
FM = "full_management"
PERIODS_FM = (1996, 2010)


def fm_log(src, start):
    sub = [] if start == PERIOD else [grid.window(start)]
    return LOGS_DIR.joinpath(*sub, group(src, FM), f"{grid.label(fm_settings())}.log")


def fm_settings():
    return {**BASE, "wave_angle_high_fraction": HEADLINE}


def fm_launch(cell):
    src, start = cell
    log = fm_log(src, start)
    if grid.finished(log):
        return
    log.parent.mkdir(parents=True, exist_ok=True)
    e = env(src, FM, fm_settings())
    e["HAT_START_YEAR"] = str(start)
    t0 = time.perf_counter()
    p = subprocess.run([sys.executable, str(grid.HINDCAST)], env=e, cwd=str(grid.PROJECT_ROOT),
                       capture_output=True, text=True, encoding="utf-8", errors="replace",
                       timeout=grid.RUN_TIMEOUT_S)
    log.write_text((p.stdout or "") + "\n--- STDERR ---\n" + (p.stderr or ""), encoding="utf-8")
    what = "done" if p.returncode == 0 else f"FAILED ({step2.stop_reason(log)})"
    print(f"{what} {group(src, FM)} {grid.window(start)} in "
          f"{(time.perf_counter() - t0) / 60:.1f} min", flush=True)


def cmd_run_2010(a):
    """The 2010-2024 full-management pair at the study's headline setting."""
    grid.check_barrier3d()
    common.keep_awake()
    todo = [(src, 2010) for src in SOURCES if not grid.finished(fm_log(src, 2010))]
    print(f"{len(todo)} to run", flush=True)
    with ThreadPoolExecutor(max_workers=a.jobs) as pool:
        list(pool.map(fm_launch, todo))
    return cmd_plot()


def fm_rates(src, start):
    """The full-management run at the headline setting, or None if not run."""
    import pandas as pd
    root = STUDY_DIR / "runs" / group(src, FM) / grid.window(start) / "zeroBE"
    for md in sorted(root.glob("*/*_run_metadata.json")):
        m = json.loads(md.read_text(encoding="utf-8"))
        if abs(float(m["wave climate"]["wave_angle_high_frac"]) - HEADLINE) < 1e-9:
            return pd.read_csv(md.parent / "tables" / "shoreline_change_rate.csv"
                               ).set_index("gis_domain")
    return None


def cmd_plot(_=None):
    """This study's figures through house_figures: (a) full management
    1996-2010, (b) full management 2010-2024, at the headline setting."""
    panels = []
    for start in PERIODS_FM:
        a, b = grid.window(start).split("_")
        panels.append(dict(label=f"Full management {a}–{b}", start=int(a), end=int(b),
                           rates={src: fm_rates(src, start) for src in SOURCES}))
    setting = (f"Hs {BASE['hs']:g} m, Tp {BASE['wave_period_s']:g} s, asymmetry "
               f"{BASE['wave_asymmetry']:g}, high-angle {HEADLINE:g}")
    note = (f"(a) full management 1996-2010, (b) full management 2010-2024; waves {setting} "
            f"({WAVE_NOTE}); island offset {OFFSET_NOTE}; no source/sink correction at the "
            f"ends (zeroBE); relocations and groins off. The natural runs are in "
            f"tables/all_runs.csv. {FORM_NOTE}")
    for p in house_figures(panels, FIG, LEGEND_TITLE, note, "full_management"):
        print(p.relative_to(STUDY_DIR))
    return 0


# Observations for the house-form figures: smoothed at 7 domains, the research
# group's range (Hannah, 2026-09-28), the southern 10 domains left raw.
LOESS_DOMAINS = 7
SKIP_SOUTHERN = 10


def coastsat_target_loess(start, window):
    """The CoastSat LRR target built as the runner builds it, at LOESS_DOMAINS,
    the rate fitted on `window` ("1996_2010", or "1996_2024" for the long-term rate)."""
    from site_layer.hat_observed_rates import COASTSAT_LRR_ROOT
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS
    from cascade_pipeline.hindcast import build_target_table
    from cascade_pipeline.coastsat_loess import (CoastSatDataset, LoessConfig,
                                                 build_coastsat_series)
    ds = CoastSatDataset(label=f"CoastSat LRR ({window.replace('_', '-')})",
                         period_start=start,
                         csv_path=str(COASTSAT_LRR_ROOT / window / "transect_lrr_full.csv"))
    cfg = LoessConfig(window_domains=(LOESS_DOMAINS,), skip_southern_domains=SKIP_SOUTHERN)
    cs = build_coastsat_series([ds], active_period_start=start, loess_config=cfg,
                               domains=HATTERAS_DOMAINS)[0]
    return build_target_table(cs, cfg, HATTERAS_DOMAINS, LOESS_DOMAINS).set_index(
        "gis_domain")["target_lrr_m_yr"]


def smooth_loess(series):
    """A per-domain series smoothed as the target is: LOESS over LOESS_DOMAINS,
    the southern SKIP_SOUTHERN left raw."""
    import pandas as pd
    from statsmodels.nonparametric.smoothers_lowess import lowess
    x = series.index.to_numpy(dtype=float)
    y = series.to_numpy(dtype=float)
    ok = np.isfinite(y)
    out = pd.Series(np.nan, index=series.index)
    out[ok] = lowess(y[ok], x[ok], frac=LOESS_DOMAINS / len(x), return_sorted=False)
    raw = series.index <= SKIP_SOUTHERN
    out[raw] = series[raw]
    return out


def duneline_change(start, end):
    """Observed dune-line net change per domain (m), between the lines that bound the window."""
    import pandas as pd
    from site_layer.hat_observed_rates import DUNELINE_ENDPOINT_ROOT
    t = pd.read_csv(DUNELINE_ENDPOINT_ROOT / f"{start}_{end}" / "domain_endpoint_summary.csv")
    return t.set_index("domain_number")["mean_change_m"]


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--jobs", type=int, default=8)
    sub.add_parser("score")
    sub.add_parser("plot")
    r = sub.add_parser("run-2010")
    r.add_argument("--jobs", type=int, default=2)
    a = ap.parse_args()
    return {"run": cmd_run, "score": cmd_score, "plot": cmd_plot,
            "run-2010": cmd_run_2010}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
