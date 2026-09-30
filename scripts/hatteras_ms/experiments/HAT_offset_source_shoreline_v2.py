"""Shoreline offset v1 vs v2, beside the dune line, on the adopted setup (2026-09-29).

Asked by Hannah on 2026-09-29: re-run the shoreline arm on shoreline offset
v2, the CoastSat mean over +/-1 yr of the start DEM's lidar flights
(1995-10-12..1997-10-12 for 1996, 2008-08-17..2010-08-17 for 2010), which
became CURRENT that day. Every shoreline run before then read v1 (the
calendar means, 1995-1997 and 2009-2011).

The option A study (2026-09-28) ran before three adopted changes -- the
per-cell dune ceilings and the beach/dune cap fix (09-28, later that day)
and the split12 storm files (09-29) -- so its v1 runs cannot be set beside
a v2 run made now. Hannah chose a clean three-way study instead: all three
offsets on today's setup, so v1 -> v2 is the only difference between the
shoreline arms.

    arms      duneline     2-brie-offset/<year>/duneline/CURRENT (v1)
              shoreline_v1 2-brie-offset/<year>/shoreline/v1, pinned by
                           HAT_OFFSET_VERSION_<year>_SHORELINE=v1
              shoreline_v2 2-brie-offset/<year>/shoreline/v2 (CURRENT)
    waves     option A (Hs 2.0 m, Tp 7.5 s, asymmetry 0.6, high-angle 0.5)
    ends      zeroBE, as in the option A study (edgeBE ends were solved on
              the dune-line offset and would favour it)
    scope     natural and full management, 1996-2010 and 2010-2024;
              relocations and groins off; 12 runs
    setup     everything else is the code default: dune ceilings, storm
              series (v3_split12_trim24), beach/dune cap

Every run's metadata is checked for the offset it actually read.

SCORES, as in the option A study
    vs CoastSat   each period's own CoastSat LRR, LOWESS 7 domains: raw and
                  smoothed share of the alongshore variation explained,
                  interior GIS 2-89, bias and r (the option A headline)
    own feature   dune-line arm against dune-line net change (m); shoreline
                  arms against total (own LRR x 14 yr) and projected
                  (1996-2024 LRR x 14 yr) shoreline change

WHERE: output/raw_runs/experiments/island-offset/2026-09-29-shoreline-offset-v1-vs-v2-adopted-setup/

    python scripts/hatteras_ms/experiments/HAT_offset_source_shoreline_v2.py run --jobs 4
    python scripts/hatteras_ms/experiments/HAT_offset_source_shoreline_v2.py score
    python scripts/hatteras_ms/experiments/HAT_offset_source_shoreline_v2.py plot
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

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_offset_source_comparison as base  # noqa: E402
import HAT_offset_source_comparison_option_a as oa  # noqa: E402

grid, common, step2 = base.grid, base.common, base.step2

TAG = "island-offset/2026-09-29-shoreline-offset-v1-vs-v2-adopted-setup"
STUDY_DIR = grid.RAW_RUNS / "experiments" / TAG
TABLES, LOGS, FIG = STUDY_DIR / "tables", STUDY_DIR / "logs", STUDY_DIR / "figures"

SETTINGS = dict(oa.OWN_SETTINGS)          # option A, high-angle 0.5
PERIODS = (1996, 2010)
SCENARIOS = ("natural", "full_management")
ARMS = {
    # arm: (HAT_ISLAND_OFFSET_SOURCE, pinned version or None, expected metadata)
    "duneline": ("duneline", None, "duneline/v1"),
    "shoreline_v1": ("shoreline", "v1", "shoreline/v1"),
    "shoreline_v2": ("shoreline", None, "shoreline/v2"),
}
WINDOWS = {"shoreline_v1": "over calendar 1995-1997 for 1996 and 2009-2011 for 2010 (v1)",
           "shoreline_v2": ("over 1995-10-12 to 1997-10-12 for 1996 and 2008-08-17 to 2010-08-17 "
                            "for 2010, +/-1 yr of the start DEM's lidar flights (v2)")}
YEARS = oa.YEARS


def cells():
    return list(product(ARMS, SCENARIOS, PERIODS))


def group(arm, scenario):
    return f"{arm}_{scenario}"


def log_path(arm, scenario, start):
    return LOGS / oa.period_label(start) / f"{group(arm, scenario)}.log"


def env(arm, scenario, start):
    src, pin, _ = ARMS[arm]
    e = grid.run_env("x", scenario, start, SETTINGS)
    e["HAT_ISLAND_OFFSET_SOURCE"] = src
    if pin:
        e[f"HAT_OFFSET_VERSION_{start}_{src.upper()}"] = pin
    e["HAT_RUN_TAG"] = f"{TAG}/runs/{group(arm, scenario)}"
    return e


def launch(cell):
    arm, scenario, start = cell
    log = log_path(*cell)
    if grid.finished(log):
        return
    log.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    p = subprocess.run([sys.executable, str(grid.HINDCAST)], env=env(*cell),
                       cwd=str(grid.PROJECT_ROOT), capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=grid.RUN_TIMEOUT_S)
    log.write_text((p.stdout or "") + "\n--- STDERR ---\n" + (p.stderr or ""), encoding="utf-8")
    what = "done" if p.returncode == 0 else f"FAILED ({step2.stop_reason(log)})"
    print(f"{what} {group(arm, scenario)} {oa.period_label(start)} in "
          f"{(time.perf_counter() - t0) / 60:.1f} min", flush=True)


def cmd_run(a):
    grid.check_barrier3d()
    common.keep_awake()
    todo = [c for c in cells() if not grid.finished(log_path(*c))]
    print(f"{len(cells())} runs, {len(todo)} to run, {a.jobs} at a time", flush=True)
    with ThreadPoolExecutor(max_workers=a.jobs) as pool:
        list(pool.map(launch, todo))
    return cmd_score()


def run_dir(arm, scenario, start):
    """The one run of this arm, checked against the offset it read."""
    root = STUDY_DIR / "runs" / group(arm, scenario) / oa.period_label(start) / grid.PRESET
    dirs = [d for d in sorted(root.glob("*")) if d.is_dir()] if root.is_dir() else []
    if len(dirs) != 1:
        raise SystemExit(f"{root}: expected one run, found {len(dirs)}")
    d = dirs[0]
    md = json.loads(next(d.glob("*_run_metadata.json")).read_text(encoding="utf-8"))
    got = str(md["identity"]["island_offset_version"])
    if got != ARMS[arm][2]:
        raise SystemExit(f"{d}: ran on offset {got!r}, expected {ARMS[arm][2]!r}")
    return d


def rates(arm, scenario, start):
    import pandas as pd
    return pd.read_csv(run_dir(arm, scenario, start) / "tables" / "shoreline_change_rate.csv"
                       ).set_index("gis_domain")


def cmd_score(_=None):
    import pandas as pd
    rows = []
    for arm, scenario, start in cells():
        rec = dict(arm=arm, scenario=scenario, period=oa.period_label(start))
        log = log_path(arm, scenario, start)
        try:
            d = run_dir(arm, scenario, start)
        except SystemExit as e:
            rec["status"] = step2.stop_reason(log) if log.is_file() else f"not run ({e})"
            rows.append(rec)
            continue
        rec["status"] = "scored"
        # vs CoastSat, the option A headline
        rec.update(grid.score_run(d, oa.coastsat_target7(start)))
        # own feature, metres over 14 yr
        rt = pd.read_csv(d / "tables" / "shoreline_change_rate.csv").set_index("gis_domain")
        if arm == "duneline":
            feats = {"own_dune_change": (oa.smooth7(oa.duneline_change(start)),
                                         rt.change_rate_m_yr * YEARS)}
        else:
            feats = {"own_total_change": (oa.coastsat_target7(start) * YEARS, rt.lrr_m_yr * YEARS),
                     "own_projected_change": (oa.coastsat_target7(1996, oa.LONG_WINDOW) * YEARS,
                                              rt.lrr_m_yr * YEARS)}
        for k, (obs, mod) in feats.items():
            s = common.alongshore_scores(mod, obs)
            mi, oi = common.interior(mod), common.interior(obs.reindex(mod.index))
            rec[f"{k}_explained"] = s["variance_explained"]
            rec[f"{k}_r"] = s["r_alongshore"]
            rec[f"{k}_bias_m"] = float((mi - oi).mean())
        rec["run_dir"] = str(d.relative_to(STUDY_DIR)).replace("\\", "/")
        rows.append(rec)
    t = pd.DataFrame(rows)
    TABLES.mkdir(parents=True, exist_ok=True)
    t.to_csv(TABLES / "all_runs.csv", index=False)
    show = ["arm", "scenario", "period", "status", "raw_variance_explained",
            "smoothed_variance_explained", "bias_m_yr", "raw_r"]
    print(t[[c for c in show if c in t]].round(3).to_string(index=False))
    diff_table()
    return 0


def diff_table():
    """How far v2 moves the model from v1: total change (LRR x 14 yr), per run pair."""
    import pandas as pd
    rows = []
    for scenario, start in product(SCENARIOS, PERIODS):
        try:
            d = (rates("shoreline_v2", scenario, start).lrr_m_yr
                 - rates("shoreline_v1", scenario, start).lrr_m_yr) * YEARS
        except SystemExit:
            continue
        di = common.interior(d)
        rows.append(dict(scenario=scenario, period=oa.period_label(start),
                         mean_m=float(di.mean()), mean_abs_m=float(di.abs().mean()),
                         max_abs_m=float(di.abs().max()), at_gis=int(di.abs().idxmax())))
    if rows:
        t = pd.DataFrame(rows)
        t.to_csv(TABLES / "shoreline_v2_minus_v1_total_change.csv", index=False)
        print("\nv2 minus v1, model total change (interior, m):")
        print(t.round(2).to_string(index=False))


def cmd_plot(_=None):
    """The house figures (full management, both periods), once per shoreline
    version, and the v2-minus-v1 difference."""
    FIG.mkdir(parents=True, exist_ok=True)
    note_base = ("Option A waves (Hs 2.0 m, Tp 7.5 s, asymmetry 0.6, high-angle 0.5); island "
                 "offset in metres; no source/sink correction at the ends (zeroBE); relocations "
                 "and groins off; adopted dune ceilings and split12 storm series.")
    for ver in ("v1", "v2"):
        panels = []
        for start in PERIODS:
            a, b = oa.period_label(start).split("_")
            panels.append(dict(label=f"Full management {a}–{b}", start=int(a), end=int(b),
                               rates={"duneline": rates("duneline", "full_management", start),
                                      "shoreline": rates(f"shoreline_{ver}", "full_management",
                                                         start)}))
        note = (f"(a) full management 1996-2010, (b) full management 2010-2024. Shoreline offset "
                f"{ver} (duneline/v1). " + note_base)
        for p in base.house_figures(panels, FIG, "No source/sink correction at the ends", note,
                                    f"full_management_shoreline_{ver}",
                                    shoreline_window=WINDOWS[f"shoreline_{ver}"]):
            print(p.relative_to(STUDY_DIR))
    print(diff_figure(note_base).relative_to(STUDY_DIR))
    return 0


def diff_figure(note_base):
    """Model total change on shoreline v2 minus v1, per domain: natural and
    full management, both periods, on the difference figures' fixed axis."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator
    from site_layer.hat_figure_style import (DOMAIN_AXIS_LABEL, INK_MUTED, C_1984, C_1997,
                                             apply_style, open_frame, record_caption, save,
                                             town_bands, _title)
    apply_style()
    fig, axes = plt.subplots(2, 2, figsize=(7.48, 4.6), sharex=True, sharey=True,
                             constrained_layout=True)
    for i, scenario in enumerate(SCENARIOS):
        for j, start in enumerate(PERIODS):
            ax = axes[i, j]
            d = (rates("shoreline_v2", scenario, start).lrr_m_yr
                 - rates("shoreline_v1", scenario, start).lrr_m_yr) * YEARS
            ax.bar(d.index, d.values, width=0.85,
                   color=[C_1997 if v > 0 else C_1984 for v in d.values])
            ax.axhline(0, color=INK_MUTED, lw=0.7)
            ax.set_xlim(0.5, 90.5)
            ax.set_ylim(-20, 30)
            ax.yaxis.set_major_locator(MultipleLocator(10))
            ax.grid(axis="y")
            open_frame(ax)
            town_bands(ax)
            a, b = oa.period_label(start).split("_")
            _title(ax, 2 * i + j, f"{scenario.replace('_', ' ').capitalize()} {a}–{b}")
            if i == 1:
                ax.set_xlabel(DOMAIN_AXIS_LABEL)
            if j == 0:
                ax.set_ylabel("v2 minus v1 (m)")
    png = FIG / "total_change_difference_shoreline_v2_minus_v1.png"
    save(fig, png, dpi=300, close=True)
    record_caption(png, (
        "Modelled total shoreline change (each run's own LRR x 14 yr) started from shoreline "
        "offset v2 (CoastSat mean over +/-1 yr of the start DEM's lidar flights) minus the "
        "same run started from v1 (calendar 1995-1997 / 2009-2011 means), per domain. Blue: v2 "
        "gives the more accretional change; red: v1 does. Rows: natural, full management; "
        "columns: 1996-2010, 2010-2024. Same axis as the shoreline-minus-dune-line difference "
        "figures (-20 to 30 m). " + note_base))
    return png


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--jobs", type=int, default=4)
    sub.add_parser("score")
    sub.add_parser("plot")
    a = ap.parse_args()
    return {"run": cmd_run, "score": cmd_score, "plot": cmd_plot}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
