"""
Do dune management heights measured from the lidar, instead of the defaults, fix the full-management dunes?

    python scripts/real-world-cases/Hatteras/dune-storm-response/HAT_measured_dune_heights.py run
    python scripts/real-world-cases/Hatteras/dune-storm-response/HAT_measured_dune_heights.py score

Three heights, all m NAVD88, from survey_dune_crests.csv (2001's unclassified survey left out):

    village dune cap      8.0  BeachDuneManager cuts every village dune cell to 4 m above the berm each year
                               (5.7 m NAVD88, obx-cascade's Nags Head default) and the model's village dunes
                               sit on it in both windows; the lidar's village crests reach 7.5-8.2 m (p95)
    road design height    5.4  the 2009->2014 rebuilds on the road domains (GIS 11-13, 51-53, 88) reached
                               4.6-6.7 m, median 5.4; was 4.3 (Velasquez-Montoya et al. 2021)
    road rebuild trigger  3.5  the same rebuilds started from 1.9-3.5 m; was the berm floor (berm + 0.3 m)

RoadwayManager resets every dune cell of a domain to the design height once any front-row cell is below
the trigger, so a higher trigger also cuts tall cells down. Two variants, both windows, managed edgeBE
with the groin as the matrix runs:

    control             the code as it is
    cap_only            the village cap alone, the road rebuild as the setup

A variant name may carry "_Hs<m>" for a wave height other than the code default (2.0 m): control_Hs3 and
cap_only_Hs3 rerun the comparison at Hs 3 m, each in its own folder.
    measured           the three heights, the stock rebuild
    measured_raise_only the three heights, the rebuild raising cells to the design height but never
                        lowering them (as HAT_road_dune_rebuild.py)

The run metadata still prints the setup's 3.94 / 0.0 m MHW rebuild heights; the EXPERIMENT line in each
log gives the ones used.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-28
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
for _path in (_HAT_DIR, _HAT_DIR / "cascade_pipeline" / "plotting"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))
from hatteras_paths import WORKSPACE as PROJECT_ROOT  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
TAG = "2026-10-28-measured-dune-heights"
EXP_DIR = PROJECT_ROOT / "output" / "raw_runs" / "dune-storm-response" / TAG
HINDCAST = _HAT_DIR / "hatteras_ms" / "HAT_hindcast.py"

WINDOWS = {1996: 2009, 2009: 2025}
PRESET, SCENARIO = "edgeBE", "full_management"
MHW_NAVD88_M = 0.36
BERM_NAVD88_M = 1.7
# m NAVD88
VILLAGE_CAP = 8.0
ROAD_DESIGN = 5.4
ROAD_TRIGGER = 3.5
VARIANTS = ("control", "cap_only", "measured", "measured_raise_only")
WORKERS = 3
# -----------------------------------------------------------------------------


# A variant name with an "_Hs<m>" suffix runs at that wave height (the code default otherwise), in its own
# folder: "cap_only_Hs3" is cap_only at Hs 3 m
def split_hs(variant):
    base, _, hs = variant.partition("_Hs")
    return base, (float(hs) if hs else None)


# Child process: the variant's heights patched in, then the unchanged hindcast
def _launch(variant):
    import runpy
    variant, _ = split_hs(variant)
    if variant != "control":
        import cascade.beach_dune_manager as BDM
        import cascade.roadway_manager as RM
        import cascade_pipeline.hindcast as H

        # The village cap, m above the berm, set on every manager after obx-cascade's own __init__
        orig_init = BDM.BeachDuneManager.__init__

        def init(self, *a, **k):
            orig_init(self, *a, **k)
            self._artificial_maximum_dune_height = VILLAGE_CAP - BERM_NAVD88_M

        BDM.BeachDuneManager.__init__ = init

        if variant == "measured_raise_only":
            orig_rebuild = RM.rebuild_dunes

            # Cells below the design come up to it; cells above keep their height
            def rebuild_dunes(yxz_dune_grid, *a, **k):
                new, _ = orig_rebuild(yxz_dune_grid, *a, **k)
                new = np.maximum(yxz_dune_grid, new)
                return new, np.sum(new - yxz_dune_grid)

            RM.rebuild_dunes = rebuild_dunes

        orig_build = H.build_cascade

        def build_cascade(*a, **k):
            n = len(k["dune_design_elevation"])
            k.update(dune_design_elevation=[ROAD_DESIGN - MHW_NAVD88_M] * n,
                     dune_minimum_elevation=[ROAD_TRIGGER - MHW_NAVD88_M] * n)
            print(f"\nEXPERIMENT {variant}: village cap {VILLAGE_CAP} m NAVD88 "
                  f"({VILLAGE_CAP - BERM_NAVD88_M:.2f} m above the berm), road design {ROAD_DESIGN} m and "
                  f"trigger {ROAD_TRIGGER} m NAVD88"
                  + (", raise-only rebuild" if variant == "measured_raise_only" else ""), flush=True)
            return orig_build(*a, **k)

        if variant != "cap_only":
            H.build_cascade = build_cascade
        else:
            print(f"\nEXPERIMENT cap_only: village cap {VILLAGE_CAP} m NAVD88 "
                  f"({VILLAGE_CAP - BERM_NAVD88_M:.2f} m above the berm), road rebuild as the setup", flush=True)
    sys.path.insert(0, str(HINDCAST.parent))
    sys.argv = [str(HINDCAST)]
    runpy.run_path(str(HINDCAST), run_name="__main__")


def run_dir(variant, start):
    base = EXP_DIR / "runs" / variant / f"{start}_{WINDOWS[start]}" / PRESET
    hits = sorted(d for d in base.glob("HAT_*") if any(d.glob("*_run_metadata.json"))) if base.exists() else []
    return hits[-1] if hits else None


# One run in its own process, logged, with its launch record
def launch(job):
    variant, start = job
    logs = EXP_DIR / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    hs = split_hs(variant)[1]
    if hs is not None:
        env["HAT_HS"] = repr(hs)
    env.update({"HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(start),
                "HAT_SOURCE_SINK_PRESET": PRESET, "HAT_SCENARIO": SCENARIO, "HAT_GROIN_ENABLED": "1",
                "HAT_RUN_KIND": "dune-storm-response", "HAT_RUN_TAG": f"{TAG}/runs/{variant}",
                "HAT_SAVE_MODEL_STATE": "1", "HAT_MAKE_GIFS": "0", "HAT_SHOW_FIGURES": "0",
                "HAT_OVERWRITE": "1", "MPLBACKEND": "Agg", "PYTHONIOENCODING": "utf-8"})
    t0 = time.time()
    log = logs / f"{variant}_{start}.log"
    with open(log, "w", encoding="utf-8") as fh:
        p = subprocess.run([sys.executable, str(_HERE), "_launch", variant],
                           env=env, cwd=HINDCAST.parent, stdout=fh, stderr=subprocess.STDOUT)
    rec = dict(variant=variant, start=start, hs=hs, village_cap=VILLAGE_CAP, road_design=ROAD_DESIGN,
               road_trigger=ROAD_TRIGGER, returncode=p.returncode, minutes=round((time.time() - t0) / 60, 1))
    with open(logs / "launches.jsonl", "a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec) + "\n")
    print(f"{variant:20s} {start}  exit {p.returncode}  {rec['minutes']} min", flush=True)


# Per run: the dune misfit against the lidar (all, villages, road), the shoreline against CoastSat, and
# how the two misfits line up, over the interior domains
def score():
    import dune_crest_comparison as dc
    from cascade_pipeline.run_layout import TABLES
    from site_layer.hatteras_site_config import HATTERAS_COMMUNITY_ZONES, SCORE_INTERIOR_GIS

    lo, hi = SCORE_INTERIOR_GIS
    village = {g for a, b in HATTERAS_COMMUNITY_ZONES for g in range(a, b + 1)}
    rows = []
    for start, end in WINDOWS.items():
        present = sorted(p.name for p in (EXP_DIR / "runs").iterdir() if p.is_dir())
        for variant in [v for v in VARIANTS if v in present] + [v for v in present if v not in VARIANTS]:
            d = run_dir(variant, start)
            if d is None:
                print(f"  {variant} {start}: no run")
                continue
            t = pd.read_csv(d / TABLES / dc.TABLE)
            t = t[(t.survey > start) & t.gis.between(lo, hi)]
            road = pd.read_csv(d / TABLES / "road_management.csv")
            on_road = set(road.gis)
            model, observed = dc.position_change(d, (start, end))
            gap = (model - observed).loc[lo:hi]
            misfit = t.groupby("gis").model_minus_survey_m.mean()
            pair = pd.DataFrame({"misfit": misfit, "gap": gap}).dropna()
            m = t.model_minus_survey_m
            rows.append(dict(
                window=f"{start}_{end}", variant=variant,
                dune_bias=m.mean(), dune_rmse=np.sqrt((m ** 2).mean()),
                village_bias=t[t.gis.isin(village)].model_minus_survey_m.mean(),
                road_bias=t[t.gis.isin(on_road)].model_minus_survey_m.mean(),
                within_half_m_pct=(m.abs() <= 0.5).mean() * 100,
                shore_bias=gap.mean(), shore_rmse=np.sqrt((gap ** 2).mean()),
                r_misfit_gap=pair.misfit.corr(pair.gap),
                road_rebuilds=int(road.dunes_rebuilt.sum()), run=d.name))
    out = pd.DataFrame(rows)
    (EXP_DIR / "tables").mkdir(parents=True, exist_ok=True)
    out.round(3).to_csv(EXP_DIR / "tables" / "scores.csv", index=False)
    print(out.drop(columns="run").round(2).to_string(index=False))


# Shoreline position change over each window, the original setup against the village cap, under CoastSat,
# with the cap's change beneath each window
def figure(suffix=""):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import dune_crest_comparison as dc
    from cascade_pipeline.run_layout import TABLES
    from site_layer.hat_figure_style import (DOMAIN_AXIS_LABEL, INK, INK_MUTED, apply_style, open_frame,
                                             record_caption, save, structures, town_bands)

    apply_style()
    fig, axes = plt.subplots(4, 1, figsize=(11, 9.5), sharex=True, layout="constrained",
                             gridspec_kw=dict(height_ratios=(3, 1.1, 3, 1.1)))
    letters = iter("abcd")
    for k, (start, end) in enumerate(WINDOWS.items()):
        ax, dax = axes[2 * k], axes[2 * k + 1]
        control, observed = dc.position_change(run_dir("control" + suffix, start), (start, end))
        capped, _ = dc.position_change(run_dir("cap_only" + suffix, start), (start, end))
        ax.axhline(0, color=INK_MUTED, lw=0.6)
        ax.plot(observed.index, observed.values, color="black", lw=2.0, label="CoastSat")
        ax.plot(control.index, control.values, color="#E69F00", lw=1.7, label="Original setup (cap 5.7 m)")
        ax.plot(capped.index, capped.values, color="#1B7F6B", lw=1.7, ls=(0, (5, 2)),
                label=f"Village cap {VILLAGE_CAP:g} m NAVD88")
        ax.set_ylabel(f"Position change,\n{start}–{end} (m)\n← landward   seaward →")
        hs_text = f", Hs {split_hs('x' + suffix)[1]:g} m" if suffix else ""
        ax.set_title(f"({next(letters)}) {start}–{end}, full management{hs_text}", loc="left", fontweight="bold")
        # The bars unsmoothed, so a village's change is not spread onto the road domains beside it
        raw = {v: pd.read_csv(run_dir(v + suffix, start) / TABLES / "shoreline_change_rate.csv")
               .set_index("gis_domain").change_rate_m_yr * (end - start) for v in ("control", "cap_only")}
        diff = raw["cap_only"] - raw["control"]
        dax.axhline(0, color=INK_MUTED, lw=0.6)
        dax.bar(diff.index, diff.values, color=np.where(diff.values >= 0, "#2166ac", "#b2182b"), width=0.8)
        dax.set_ylabel("Cap − original\n(m)")
        dax.set_title(f"({next(letters)}) change from the village cap, {start}–{end}", loc="left",
                      fontweight="bold")
        lim = max(1.0, float(np.nanmax(np.abs(diff.values))) * 1.15)
        dax.set_ylim(-lim, lim)
        for a in (ax, dax):
            open_frame(a)
            a.set_xlim(0.5, 90.5)
        town_bands(ax)
        town_bands(dax, label=False)
        structures(ax, label=k == 0)
    axes[0].legend(loc="lower center", frameon=False, ncol=3)
    axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)
    png = EXP_DIR / "figures" / f"position_change_village_cap_vs_original{suffix}.png"
    save(fig, png, dpi=300, close=True)
    record_caption(png, (
        "Shoreline position change over each window under full management (edgeBE, groin on): CoastSat "
        "(black), the original setup (orange; BeachDuneManager caps village dunes at 4 m above the berm, "
        f"5.7 m NAVD88) and the village cap raised to {VILLAGE_CAP:g} m NAVD88 from the lidar (green, "
        "dashed; the road rebuild as the original). Both model lines and CoastSat are 7-domain LOWESS, as on "
        "the v2 dune crest overview; seaward positive. Bars: the capped run minus the original, unsmoothed "
        "(endpoint rate x the window's years), blue seaward. Grey bands: the villages, where the beach and "
        "dune manager, and so the cap, acts. Runs: dune-storm-response/" + TAG + f"/runs/{{control,cap_only}}{suffix}"
        + (f", wave height Hs {split_hs('x' + suffix)[1]:g} m." if suffix else ", wave height the code default.")))
    print(f"  -> {png}")


def main():
    if len(sys.argv) > 2 and sys.argv[1] == "_launch":
        return _launch(sys.argv[2])
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("action", choices=("run", "score", "figure"))
    ap.add_argument("--variants", default=",".join(VARIANTS))
    ap.add_argument("--hs", default="", help='figure: the Hs suffix of the runs to draw, e.g. "_Hs3"')
    a = ap.parse_args()
    if a.action == "run":
        jobs = [(v, s) for s in WINDOWS for v in a.variants.split(",")]
        with ThreadPoolExecutor(WORKERS) as pool:
            list(pool.map(launch, jobs))
    if a.action == "figure":
        return figure(a.hs)
    score()


if __name__ == "__main__":
    main()
