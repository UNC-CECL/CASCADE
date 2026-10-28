#!/usr/bin/env python3
"""
A tapered fillet: the groin's updrift gain spread over GIS 6, 7 and 8 (1/2, 1/3, 1/6), scored in both hindcast
windows against the runner's groin (all of it in GIS 6) and no groin, all on today's runner.

    python tapered_fillet.py run [--streams 4]   # no groin + both layouts over the strength grid, both windows
    python tapered_fillet.py score               # -> tapered_fillet_scores.csv, tapered_fillet_domains.csv

The taper is the advisor's suggestion (2026-10-20): the full effect against the groin, two-thirds of it one domain
north, one-third two domains north, scaled to sum to one so the groin still only moves sand. The downdrift loss
stays in GIS 5. Patched in-process as ../2026-10-14-face-6-7/field_windows.py; `cascade/groin.py` and the runner
are not edited. Everything else at the runner's defaults (edgeBE, full management, relocations off, off from 2017).
Scored as the field-windows study: the gap across the groins (GIS 5 minus GIS 6) against annual CoastSat, and net
change against the CoastSat net-change target over GIS 2-89 and GIS 2-12, plus net change per domain at GIS 4-9.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-20
"""

from __future__ import annotations

# Non-code files beside this script live in the workspace on D: (hatteras_paths.workspace_twin)
import sys as _hat_sys  # noqa: E402
from pathlib import Path as _HatPath  # noqa: E402
_HAT_DIR = next(_p for _p in _HatPath(__file__).resolve().parents if (_p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import workspace_twin  # noqa: E402
import argparse
import os
import queue
import subprocess
import sys
import threading
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
# Hatteras/ holds the scripts; data/ and output/ live in the workspace on D: (hatteras_paths)
import sys as _hat_sys  # noqa: E402
_HAT_DIR = next(p for p in HERE.parents  if (p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE as REPO  # noqa: E402
for p in (_HAT_DIR, _HAT_DIR / "hatteras_ms", _HAT_DIR / "hatteras_ms" / "experiments",
          _HAT_DIR / "analyze_output" / "compare_runs" / "matrix_vs_observed",
          HERE.parent / "2026-10-09-corrected-groin"):
    sys.path.insert(0, str(p))

# --- CONFIG ------------------------------------------------------------------
RUNNER = _HAT_DIR / "hatteras_ms" / "HAT_hindcast.py"
STUDY = "groin-wave-driven/2026-10-20-tapered-fillet"
RAW = REPO / "output" / "raw_runs" / "experiments" / STUDY
LOGS = workspace_twin(HERE) / "logs"
WINDOW = {1996: "1996_2009", 2009: "2009_2025"}
STRENGTHS = {1996: (0.10, 0.13, 0.16, 0.19, 0.22, 0.25), 2009: (0.06, 0.08, 0.10, 0.12, 0.16, 0.20)}
TAPER = (1 / 2, 1 / 3, 1 / 6)        # share of the updrift gain in GIS 6, 7, 8
LAYOUTS = ("face", "taper")
NET_KEY = {1996: (1996, 2009), 2009: (2009, 2025)}
STORMS = "v3_240"                     # checked in each run's log
# Dune growth bounds, pinned: the runner's RMIN/RMAX changed three times while this study ran (2026-10-20), and
# every run must share them. 0.005-0.01 is what most of the runs had; Dmaxel 12.2 m and slope 0.06 come unchanged
# from the runner. Checked against each run's parameter yaml at scoring
DUNE_GROWTH = (0.005, 0.01)
# -----------------------------------------------------------------------------


# The taper, patched into the GroinManager Cascade builds before the runner runs. Each year the gain the groin put in GIS 6
# is re-shared over GIS 6-8, in the step and in the groin's own contribution the footprint is built from
def patch_groin():
    import cascade.cascade as cc
    import cascade.groin as cg
    base = cg.GroinManager

    class TaperedGroin(base):
        def update(self, brie, x_s_dt):
            x_s_dt = super().update(brie, x_s_dt)
            gain = self.record[-1]["applied_dx_updrift_m"]
            self.record[-1]["taper_moved_m"] = 0.0
            if gain >= 0.0:                              # drift from the south this year: the gain is not in GIS 6
                return x_s_dt
            last = self._last_year
            for k, w in enumerate(TAPER[1:], start=1):   # GIS increases with pad: GIS 7 = hi + 1, GIS 8 = hi + 2
                for arr in (x_s_dt, last.dx_m, last.groin_dx_m):
                    arr[self._hi_pad] -= w * gain
                    arr[self._hi_pad + k] += w * gain
            self.record[-1]["taper_moved_m"] = -(1 - TAPER[0]) * gain
            return x_s_dt

    cc.GroinManager = cg.GroinManager = TaperedGroin


def member(layout, s):
    return "nogroin" if layout is None else f"{layout}/s{s:.2f}"


def run_dir(period, layout, s):
    hits = sorted((RAW / member(layout, s) / WINDOW[period]).glob("*/*/*_shoreline_matrix.npy"))
    return hits[0].parent if hits else None


# The Barrier3D settings a run used, from its copy of the parameter yaml: dune ceiling, beach slope, growth rates
def run_settings(rd):
    import yaml
    doc = yaml.safe_load(next(rd.glob("*-parameters.yaml")).read_text())
    return tuple((k, doc[k]) for k in ("Dmaxel", "beta", "rmin", "rmax"))


def _launch():
    import runpy
    import cascade_pipeline.hindcast as hc
    build = hc.build_cascade

    def pinned_build(*a, **kw):
        n = len(kw["rmin"])
        kw["rmin"], kw["rmax"] = [DUNE_GROWTH[0]] * n, [DUNE_GROWTH[1]] * n
        return build(*a, **kw)

    hc.build_cascade = pinned_build
    if os.environ.get("TAPER") == "1":
        patch_groin()
    sys.path.insert(0, str(RUNNER.parent))
    sys.argv = [str(RUNNER)]
    runpy.run_path(str(RUNNER), run_name="__main__")


def run(period, layout, s):
    env = {k: v for k, v in os.environ.items() if not k.startswith(("HAT_", "TAPER"))}
    env.update({
        "HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(period),
        "HAT_SOURCE_SINK_PRESET": "edgeBE", "HAT_SCENARIO": "full_management",
        "HAT_GROIN_ENABLED": "false" if layout is None else "true",
        "HAT_RUN_KIND": "experiment", "HAT_RUN_TAG": f"{STUDY}/{member(layout, s)}",
        "HAT_SAVE_MODEL_STATE": "false", "HAT_SHOW_FIGURES": "false", "HAT_MAKE_GIFS": "false",
        "HAT_OVERWRITE": "false", "MPLBACKEND": "Agg", "PYTHONIOENCODING": "utf-8",
        "TAPER": "1" if layout == "taper" else "0",
    })
    if layout is not None:
        env["HAT_GROIN_STRENGTH"] = str(s)
    LOGS.mkdir(parents=True, exist_ok=True)
    log = workspace_twin(LOGS) / f"{WINDOW[period]}_{member(layout, s).replace('/', '_')}.log"
    t0 = time.time()
    with open(log, "w", encoding="utf-8") as fh:
        code = subprocess.call([sys.executable, str(Path(__file__)), "_launch"], cwd=REPO, env=env,
                               stdout=fh, stderr=subprocess.STDOUT)
    d = run_dir(period, layout, s)
    if d is not None and layout == "taper":
        (d / "tapered_fillet.txt").write_text(
            f"patched in-process by tapered_fillet.py: the groin's updrift gain shared over GIS 6, 7, 8 as "
            f"{', '.join(f'{w:.3f}' for w in TAPER)}, at strength {s:g}, off from 2017. The runner's report "
            f"shows the GIS 5|6 groin.")
    print(f"{time.strftime('%H:%M:%S')}  {WINDOW[period]} {member(layout, s)}: exit {code} "
          f"({(time.time() - t0) / 60:.1f} min)", flush=True)
    return code


def cells():
    out = [(p, None, None) for p in WINDOW]
    out += [(p, lay, s) for p, ss in STRENGTHS.items() for lay in LAYOUTS for s in ss]
    return out


def cmd_run(a):
    todo = [c for c in cells() if run_dir(*c) is None]
    print(f"{len(cells()) - len(todo)} of {len(cells())} done; running {len(todo)} in {a.streams} streams", flush=True)
    jobs, failed = queue.Queue(), []
    for c in todo:
        jobs.put(c)

    def worker():
        while True:
            try:
                c = jobs.get_nowait()
            except queue.Empty:
                return
            if run(*c):
                failed.append(c)

    threads = [threading.Thread(target=worker) for _ in range(a.streams)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    print(f"DONE: {len(todo) - len(failed)}/{len(todo)} clean; failed {failed}", flush=True)


def cmd_score(_a):
    import corrected_groin as cgm
    import site_layer.hat_scoring as common
    from matrix_vs_observed import model_end_change
    from site_layer.hat_observed_rates import net_change_domain_csv
    rows, doms, settings = [], [], set()
    for period in WINDOW:
        gap_obs = cgm.observed_change(period)
        d = pd.read_csv(net_change_domain_csv(*NET_KEY[period])).set_index("gis_domain")
        net_obs = d["net_change_lowess7_m"]
        doms.append(dict(period=period, layout="observed", strength=None,
                         **{f"gis{g}": d.loc[g, "net_change_m"] for g in range(4, 10)}))
        for c in [c for c in cells() if c[0] == period]:
            rd = run_dir(*c)
            if rd is None:
                continue
            log = workspace_twin(LOGS) / f"{WINDOW[period]}_{member(c[1], c[2]).replace('/', '_')}.log"
            if STORMS not in log.read_text(encoding="utf-8", errors="replace"):
                raise RuntimeError(f"{log.name} did not run on the {STORMS} storms")
            settings.add(run_settings(rd))
            layout = {None: "no groin", "face": "runner (all in GIS 6)", "taper": "taper 1/2, 1/3, 1/6"}[c[1]]
            path = next(rd.glob("*_shoreline_matrix.npy"))
            raw = model_end_change(rd)
            net = common.smooth_like_target(raw)
            row = dict(period=period, layout=layout, strength=c[2], **cgm.score_annual(path, period, gap_obs))
            row.update(cgm.split_scores(path, period, gap_obs))
            for name, (lo, hi) in (("net_rmse_2_89", (2, 89)), ("net_rmse_2_12", (2, 12))):
                r = (net - net_obs).loc[lo:hi].dropna()
                row[name] = float(np.sqrt(np.mean(r ** 2)))
            rows.append(row)
            doms.append(dict(period=period, layout=layout, strength=c[2],
                             **{f"gis{g}": raw.loc[g] for g in range(4, 10)}))
    if len(settings) > 1:
        raise RuntimeError(f"the runs used different Barrier3D settings {sorted(settings)}; the runner's settings "
                           "changed during the study, so rerun the stale cells")
    t = pd.DataFrame(rows)
    t.to_csv(workspace_twin(HERE) / "tapered_fillet_scores.csv", index=False, float_format="%.3f")
    pd.DataFrame(doms).to_csv(workspace_twin(HERE) / "tapered_fillet_domains.csv", index=False, float_format="%.2f")
    pd.set_option("display.width", 220)
    cols = ["period", "layout", "strength", "rmse_annual", "bias_annual", "rmse_2010_2016", "rmse_2017_2025",
            "net_rmse_2_89", "net_rmse_2_12"]
    print(t[[c for c in cols if c in t]].round(2).to_string(index=False))
    print(pd.DataFrame(doms).round(1).to_string(index=False))


# The fillet near the groin at each layout's best gap strength: end-of-window change against the observed (top),
# and the groin's own effect, each layout less no groin (bottom); calibration window left, test right
FIG_STRENGTH = {1996: {"face": 0.16, "taper": 0.19}, 2009: {"face": 0.12, "taper": 0.12}}
FIG_GIS = np.arange(1, 15)


def cmd_figure(_a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from site_layer import hat_figure_style as fs
    from matrix_vs_observed import model_end_change
    from site_layer.hat_observed_rates import net_change_domain_csv
    styles = {None: dict(color=fs.C["BASE"], lw=1.3), "face": dict(color=fs.C["ACCENT"], lw=1.6),
              "taper": dict(color=fs.C["ADDED"], lw=1.6)}

    fs.apply_style()
    fig, axes = plt.subplots(2, 2, figsize=fs.figsize("double", height=5.6), sharex=True, sharey="row",
                             constrained_layout=True)
    for col, period in enumerate(WINDOW):
        d = pd.read_csv(net_change_domain_csv(*NET_KEY[period])).set_index("gis_domain")["net_change_m"]
        base = model_end_change(run_dir(period, None, None)).loc[FIG_GIS]
        top, bottom = axes[0, col], axes[1, col]
        for ax in (top, bottom):
            ax.axhline(0, color=fs.INK_MUTED, lw=0.6)
            ax.axvline(5.5, color=fs.C["GROIN"], lw=1.0, ls="--", zorder=1)
        top.text(5.6, 0.97, "groin", transform=top.get_xaxis_transform(), ha="left", va="top", fontsize=6.5,
                 color=fs.C["GROIN"])
        top.text(0.01, 0.03, f"{period}-{NET_KEY[period][1]}" + (" (calibration)" if period == 1996 else " (test)"),
                 transform=top.transAxes, ha="left", va="bottom", fontsize=7.5)
        # Where the taper puts the updrift gain
        bottom.axvspan(5.5, 8.5, color="0.94", zorder=0)
        top.plot(FIG_GIS, d.loc[FIG_GIS].values, "o-", ms=3, color="black", lw=1.8, zorder=5)
        top.plot(FIG_GIS, base.values, "o-", ms=2.5, **styles[None])
        for lay in LAYOUTS:
            s = FIG_STRENGTH[period][lay]
            m = model_end_change(run_dir(period, lay, s)).loc[FIG_GIS]
            top.plot(FIG_GIS, m.values, "o-", ms=2.5, **styles[lay])
            bottom.plot(FIG_GIS, (m - base).values, "o-", ms=2.5, **styles[lay])
        bottom.text(0.99, 0.97, f"strength: runner {FIG_STRENGTH[period]['face']:.2f}, "
                    f"taper {FIG_STRENGTH[period]['taper']:.2f}", transform=bottom.transAxes, ha="right",
                    va="top", fontsize=6.5, color=fs.INK_MUTED)
        bottom.set_xlabel(fs.DOMAIN_AXIS_LABEL)
        bottom.set_xticks(FIG_GIS)
    axes[0, 0].set_ylabel(f"Shoreline change over the window (m)\nseaward +")
    axes[1, 0].set_ylabel("Groin's effect: model less\nno-groin model (m), seaward +")
    axes[1, 0].text(7.0, 0.03, "taper's GIS 6-8", transform=axes[1, 0].get_xaxis_transform(), ha="center",
                    va="bottom", fontsize=6.3, color=fs.INK_MUTED)
    handles = [Line2D([], [], color="black", lw=1.8, marker="o", ms=3, label="observed (CoastSat)")]
    handles += [Line2D([], [], marker="o", ms=2.5, label=lab, **styles[lay]) for lay, lab in
                ((None, "model, no groin"), ("face", "model, runner's groin (all of the gain in GIS 6)"),
                 ("taper", "model, taper (1/2, 1/3, 1/6 of the gain in GIS 6, 7, 8)"))]
    fig.legend(handles=handles, loc="outside lower center", ncol=2, frameon=False, fontsize=6.8)
    out = fs.save(fig, workspace_twin(HERE) / "figures" / "tapered_fillet_1996_2009_2009_2025", dpi=300, close=True)
    print(*out, sep="\n")


# Plan view at 1:1, north up, GIS 4-8: the start shoreline (the observed start each run begins from), the observed
# end (start plus CoastSat's net change) and each run's end shoreline. Each domain is drawn as the flat 500 m
# segment the model treats it as. Cross-shore is the model's frame (seaward +, east), not UTM
PLAN_GIS = np.arange(4, 9)
DOMAIN_KM = 0.5


def cmd_positions(_a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from site_layer import hat_figure_style as fs
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as D
    from site_layer.hat_observed_rates import net_change_domain_csv
    pads = [D.gis_to_pad(g) for g in PLAN_GIS]
    y0 = (PLAN_GIS - 0.5 - PLAN_GIS[0] + 0.5) * DOMAIN_KM * 1000      # south edge of each domain, m
    # Widest underneath, so a run that ends where another does still shows as a rim around it
    styles = {"start": dict(color=fs.INK, lw=1.2, ls=":", zorder=5), None: dict(color=fs.C["BASE"], lw=3.4, zorder=2),
              "taper": dict(color=fs.C["ADDED"], lw=2.2, zorder=3), "face": dict(color=fs.C["ACCENT"], lw=1.1, zorder=4)}

    def seaward(period, layout, s, row):
        m = np.load(next(run_dir(period, layout, s).glob("*_shoreline_matrix.npy")))
        return -m[row, pads]

    def steps(ax, x, **kw):
        ys = np.r_[y0, y0[-1] + DOMAIN_KM * 1000]
        ax.stairs(x, ys, orientation="horizontal", baseline=None, **kw)

    fs.apply_style()
    fig, axes = plt.subplots(1, 2, figsize=(5.6, 6.6), constrained_layout=True)
    for ax, period in zip(axes, WINDOW):
        start = seaward(period, None, None, 0)
        d = pd.read_csv(net_change_domain_csv(*NET_KEY[period])).set_index("gis_domain")["net_change_m"]
        obs_end = start + d.loc[PLAN_GIS].to_numpy()
        lo = min(start.min(), obs_end.min()) - 120
        ax.fill_betweenx(np.r_[y0, y0[-1] + 500], lo, np.r_[start, start[-1]], step="post", color="0.93", lw=0)
        steps(ax, start, **styles["start"])
        steps(ax, seaward(period, None, None, -1), **styles[None])
        for lay in LAYOUTS:
            steps(ax, seaward(period, lay, FIG_STRENGTH[period][lay], -1), **styles[lay])
        ax.plot(obs_end, y0 + 250, "o", ms=5, mfc="white", mec="black", mew=1.2, zorder=6)
        # The groin, at the GIS 5|6 boundary, drawn across the start shoreline
        yg = y0[list(PLAN_GIS).index(6)]
        xg = start[list(PLAN_GIS).index(6)]
        ax.plot([xg - 150, xg + 60], [yg, yg], color=fs.C["GROIN"], lw=2.0, solid_capstyle="butt", zorder=1)
        ax.text(xg + 70, yg, "groin", color=fs.C["GROIN"], fontsize=6.5, va="center")
        for g, y in zip(PLAN_GIS, y0):
            ax.text(lo + 15, y + 250, f"GIS {g}", fontsize=6.5, color=fs.INK_MUTED, va="center")
        ax.set_aspect("equal")
        ax.set_xlim(lo, max(start.max(), obs_end.max()) + 120)
        ax.set_ylim(y0[0], y0[-1] + 500)
        ax.set_xlabel("Cross-shore position (m), seaward (east) +")
        ax.text(0.03, 0.99, f"{period}-{NET_KEY[period][1]}" + (" (calibration)" if period == 1996 else " (test)")
                + f"\nstrength: runner {FIG_STRENGTH[period]['face']:.2f},\ntaper {FIG_STRENGTH[period]['taper']:.2f}",
                transform=ax.transAxes, ha="left", va="top", fontsize=6.8)
    axes[0].set_ylabel("Alongshore distance (m), north up")
    handles = [Line2D([], [], label=f"Observed: start shoreline ({'/'.join(str(p) for p in WINDOW)}), "
                      "every model run begins here", **styles["start"]),
               Line2D([], [], ls="", marker="o", ms=5, mfc="white", mec="black", mew=1.2,
                      label="Observed: end shoreline (start + CoastSat net change)"),
               Line2D([], [], label="Model, end: no groin", **styles[None]),
               Line2D([], [], label="Model, end: runner's groin (all of the gain in GIS 6)", **styles["face"]),
               Line2D([], [], label="Model, end: taper (1/2, 1/3, 1/6 of the gain in GIS 6, 7, 8)", **styles["taper"])]
    fig.legend(handles=handles, loc="outside lower center", ncol=1, frameon=False, fontsize=6.6)
    out = fs.save(fig, workspace_twin(HERE) / "figures" / "tapered_fillet_plan_1996_2009_2009_2025", dpi=300, close=True)
    print(*out, sep="\n")


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "_launch":
        return _launch()
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sp = sub.add_parser("run")
    sp.add_argument("--streams", type=int, default=4)
    sub.add_parser("score")
    sub.add_parser("figure")
    sub.add_parser("positions")
    a = ap.parse_args()
    {"run": cmd_run, "score": cmd_score, "figure": cmd_figure, "positions": cmd_positions}[a.cmd](a)


if __name__ == "__main__":
    main()
