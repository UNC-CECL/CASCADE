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

WHERE: output/raw_runs/experiments/island-offset/2026-09-25-offset-source-duneline-vs-shoreline/
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
TAG = "island-offset/2026-09-25-offset-source-duneline-vs-shoreline"
STUDY_DIR = grid.RAW_RUNS / "experiments" / TAG
TABLES_DIR, LOGS_DIR, FIG = STUDY_DIR / "tables", STUDY_DIR / "logs", STUDY_DIR / "figures"
PERIOD = 1996
SOURCES = ("duneline", "shoreline")
SCENARIOS = grid.SCENARIOS
BASE = {"hs": 1.0, "wave_period_s": 8.0, "wave_asymmetry": 0.8}
HIGH_ANGLE = (0.3, 0.4, 0.45, 0.5, 0.55)
HEADLINE = 0.45


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
            if not str(got).startswith(src + "/"):
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


def cmd_plot(_=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd
    from matplotlib.lines import Line2D
    import HAT_metres_2_wave_sensitivity_plot as p2
    from site_layer.hat_figure_style import (INK, INK_MUTED, DOMAIN_AXIS_LABEL, _title,
                                             apply_style, open_frame, record_caption, save,
                                             structures, support_dir, town_bands)
    apply_style()
    t = pd.read_csv(TABLES_DIR / "all_runs.csv")
    target, obs = common.coastsat_target(PERIOD), p2.observed_change(PERIOD)
    col = {"duneline": "#1b7f6b", "shoreline": "#6a3d9a"}
    name = {"duneline": "Dune line", "shoreline": "Shoreline (CoastSat)"}
    out = []
    with plt.rc_context(p2.SCREEN_RC):
        # 1. the profiles at the headline setting, and their difference
        f, axes = plt.subplots(3, 2, figsize=(16, 13), sharex=True, constrained_layout=True,
                               gridspec_kw={"height_ratios": [3, 3, 2]})
        rows = []
        for j, sc in enumerate(SCENARIOS):
            prof = {}
            for src in SOURCES:
                r = t[(t.source == src) & (t.scenario == sc) & (t.status == "scored")
                      & np.isclose(t.wave_angle_high_fraction, HEADLINE)]
                if r.empty:
                    continue
                r = r.iloc[0]
                rt = pd.read_csv(STUDY_DIR / r.run_dir / "tables" / "shoreline_change_rate.csv"
                                 ).set_index("gis_domain")
                prof[src] = (rt.lrr_m_yr, rt.change_rate_m_yr * 14, r)
            for i, (o, k) in enumerate(((target, 0), (obs, 1))):
                ax = axes[i, j]
                ax.plot(o.index, o.values, color=INK, lw=2.6, zorder=6)
                for src, v in prof.items():
                    ax.plot(v[k].index, v[k].values, color=col[src], lw=1.3, alpha=0.8, zorder=4)
                    sm = common.smooth_like_target(v[k])
                    ax.plot(sm.index, sm.values, color=col[src], lw=2.6, ls=(0, (5, 3)),
                            alpha=0.6, zorder=5)
            ax = axes[2, j]
            if len(prof) == 2:
                d = prof["shoreline"][0] - prof["duneline"][0]
                ax.bar(d.index, d.values, color=[col["shoreline"] if v > 0 else col["duneline"]
                                                  for v in d.values], width=0.85)
                rows.append(dict(scenario=sc, mean_abs_diff_m_yr=float(d.abs().mean()),
                                 max_abs_diff_m_yr=float(d.abs().max()),
                                 at_gis=int(d.abs().idxmax())))
            for ax in axes[:, j]:
                ax.axhline(0, color=INK_MUTED, lw=0.6)
                ax.set_xlim(1, 90)
                ax.grid(axis="y")
                open_frame(ax)
                town_bands(ax, label=(ax is axes[0, j]), fontsize=11)
            structures(axes[1, j], label=True, label_pt=11)
            head = [f"{p2.BEST_NAME[sc]}, 1996–2010, high-angle {HEADLINE:g}"]
            for src, v in prof.items():
                r = v[2]
                head.append(f"{name[src]}: {100 * r.smoothed_variance_explained:+.0f}% "
                            f"(raw {100 * r.raw_variance_explained:+.0f}%), bias {r.bias_m_yr:+.2f}")
            _title(axes[0, j], j, "\n".join(head))
            _title(axes[1, j], 2 + j, "")
            _title(axes[2, j], 4 + j, "")
            axes[2, j].set_xlabel(DOMAIN_AXIS_LABEL)
        axes[0, 0].set_ylabel("Shoreline change rate,\nLRR (m/yr)")
        axes[1, 0].set_ylabel("Position change,\n2010 minus 1996 (m)")
        axes[2, 0].set_ylabel("Rate difference,\nshoreline − dune line (m/yr)")
        handles = [Line2D([], [], color=INK, lw=2.6, label="CoastSat (LOESS, 10 domains)")]
        for src in SOURCES:
            handles += [Line2D([], [], color=col[src], lw=1.3, label=f"Model, {name[src].lower()} offset"),
                        Line2D([], [], color=col[src], lw=2.6, ls=(0, (5, 3)), alpha=0.6,
                               label=f"{name[src]}, smoothed (scored)")]
        f.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False,
                 title="Hs 1.0 m, Tp 8 s, asymmetry 0.8, high-angle 0.45")
        png = FIG / "profiles_duneline_vs_shoreline_1996_2010.png"
        save(f, png, dpi=300, close=True)
        pd.DataFrame(rows).to_csv(support_dir(png.parent) / f"{png.stem}.csv", index=False)
        record_caption(png, (
            "The island offset (planform orientation) set from the dune line (green) and from "
            "the CoastSat shoreline (purple), 1996-2010, Hs 1.0 m, Tp 8 s, asymmetry 0.8, "
            "high-angle 0.45 (the best managed setting found, tuned with the dune-line offset). "
            "Left natural, right full management. Top: LRR rate against the CoastSat target "
            "(black); middle: position change against the observed CoastSat change; solid thin "
            "lines per domain, dashed the model smoothed like the target (the scored series). "
            "Bottom: the rate difference, shoreline minus dune line. Header: smoothed share of "
            "the alongshore variation explained (raw in brackets) and bias, m/yr."))
        out.append(png)
        # 2. score against the high-angle fraction, both offsets
        f, axes = plt.subplots(1, 3, figsize=(16, 5.6), constrained_layout=True)
        ls = {SCENARIOS[0]: "-", SCENARIOS[1]: "--"}
        for src, sc in product(SOURCES, SCENARIOS):
            x = t[(t.source == src) & (t.scenario == sc) & (t.status == "scored")
                  ].sort_values("wave_angle_high_fraction")
            kw = dict(color=col[src], ls=ls[sc], marker="o", lw=2,
                      label=f"{name[src]}, {p2.BEST_NAME[sc].lower()}")
            axes[0].plot(x.wave_angle_high_fraction, 100 * x.smoothed_variance_explained, **kw)
            axes[1].plot(x.wave_angle_high_fraction, x.bias_m_yr, **kw)
            axes[2].plot(x.wave_angle_high_fraction, x.smoothed_r, **kw)
        for ax, lab in zip(axes, ("Smoothed variation explained (%)", "Mean bias (m/yr)",
                                  "Correlation with CoastSat (smoothed)")):
            ax.axhline(0, color=INK_MUTED, lw=0.6)
            ax.set_xlabel("Fraction of high-angle waves (> 45°)")
            ax.set_ylabel(lab)
            ax.grid(True)
            open_frame(ax)
        for k, ax in enumerate(axes):
            _title(ax, k, "")
        axes[0].legend(frameon=False)
        png = FIG / "scores_vs_high_angle_duneline_vs_shoreline_1996_2010.png"
        save(f, png, dpi=300, close=True)
        t.to_csv(support_dir(png.parent) / f"{png.stem}.csv", index=False)
        record_caption(png, (
            "Scores against the high-angle fraction for the two island offsets, 1996-2010, "
            "Hs 1.0 m, Tp 8 s, asymmetry 0.8; solid natural, dashed full management. (a) "
            "share of the alongshore variation explained by the smoothed model; (b) mean bias; "
            "(c) correlation of the smoothed model with the CoastSat target. Interior GIS 2-89."))
        out.append(png)
    for p in out:
        print(p.relative_to(STUDY_DIR))
    return 0


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--jobs", type=int, default=8)
    sub.add_parser("score")
    sub.add_parser("plot")
    a = ap.parse_args()
    return {"run": cmd_run, "score": cmd_score, "plot": cmd_plot}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
