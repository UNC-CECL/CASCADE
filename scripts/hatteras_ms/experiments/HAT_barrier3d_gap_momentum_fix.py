"""
How much do the three Barrier3D overwash fixes move the results?

    python scripts/hatteras_ms/experiments/HAT_barrier3d_gap_momentum_fix.py run --workers 4
    python scripts/hatteras_ms/experiments/HAT_barrier3d_gap_momentum_fix.py compare

Four runs on the fix-branch worktree, each against its matrix twin: net
change, LRR skill and overwash against the imagery. Details: scripts/hatteras_ms/experiments/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-29
"""
from __future__ import annotations

import argparse
import importlib.util
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
PROJECT_ROOT = next(_p for _p in _HERE.parents if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

# --- CONFIG ------------------------------------------------------------------
HINDCAST = PROJECT_ROOT / "scripts" / "hatteras_ms" / "HAT_hindcast_1984_2024.py"
WORKTREE = PROJECT_ROOT.parent / "Barrier3D-overwashfix"
FIX_BRANCH = "fix/overwash-gaps-momentum"
TAG = "code-checks/2026-09-28-barrier3d-overwash-gap-momentum-fix"
EXP_DIR = PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
MATRIX = PROJECT_ROOT / "output" / "raw_runs" / "archive" / "2026-09-28-loess10-ends" / "matrix"   # the controls, archived when the matrix moved on

# member -> (start year, scenario, control run name)
MEMBERS = {
    "natural_1996": (1996, "natural", "HAT_1996_2010_edgeBE_offsetmetres_noroad_nobdm_nogroin"),
    "managed_1996": (1996, "full_management", "HAT_1996_2010_edgeBE_offsetmetres_road_bdm_nogroin"),
    "natural_2010": (2010, "natural", "HAT_2010_2024_edgeBE_offsetmetres_noroad_nobdm_nogroin"),
    "managed_2010": (2010, "full_management", "HAT_2010_2024_edgeBE_offsetmetres_road_bdm_nourish_nogroin"),
}
# -----------------------------------------------------------------------------


# A start year's window label
def window(start):
    return f"{start}_{2010 if start == 1996 else 2024}"


# A member's matrix run
def control_dir(member):
    start, _, name = MEMBERS[member]
    return MATRIX / window(start) / "edgeBE" / name


# A member's run on the fix branch
def fixed_dir(member):
    start, _, name = MEMBERS[member]
    return EXP_DIR / "runs" / f"fixed_{member}" / window(start) / "edgeBE" / name


# The environment for a member's fixed run
def env_for(member):
    start, scenario, _ = MEMBERS[member]
    env = dict(os.environ)
    env.update({
        "HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(start),
        "HAT_SOURCE_SINK_PRESET": "edgeBE", "HAT_SCENARIO": scenario,
        "HAT_RUN_KIND": "experiment", "HAT_RUN_TAG": f"{TAG}/runs/fixed_{member}",
        "HAT_SAVE_MODEL_STATE": "1", "HAT_MAKE_GIFS": "0", "HAT_SHOW_FIGURES": "0",
        "PYTHONPATH": str(WORKTREE) + os.pathsep + env.get("PYTHONPATH", ""),
        "PYTHONIOENCODING": "utf-8",
    })
    return env


# One fixed run, logged, refusing any other Barrier3D
def run_member(member):
    logs = EXP_DIR / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    log = logs / f"fixed_{member}.log"
    t0 = time.time()
    with open(log, "w", encoding="utf-8") as fh:
        p = subprocess.run([sys.executable, str(HINDCAST)], env=env_for(member),
                           cwd=HINDCAST.parent, stdout=fh, stderr=subprocess.STDOUT)
    text = log.read_text(encoding="utf-8", errors="replace")
    b3d = next((ln for ln in text.splitlines() if ln.startswith("BARRIER3D =")), "")
    ok_branch = FIX_BRANCH in b3d
    rec = dict(member=member, returncode=p.returncode, minutes=round((time.time() - t0) / 60, 1),
               barrier3d=b3d.strip(), on_fix_branch=ok_branch, run_dir=str(fixed_dir(member)))
    with open(logs / "launches.jsonl", "a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec) + "\n")
    if not ok_branch:
        raise SystemExit(f"{member}: ran on {b3d!r}, not {FIX_BRANCH}; see {log}")
    print(f"fixed {member}: exit {p.returncode}, {rec['minutes']} min, {b3d}")
    return rec


# Comparison

# overwash_vs_model.py, loaded from its folder
def _overwash_module():
    path = PROJECT_ROOT / "scripts" / "input_prep" / "8-overwash-analysis" / "4-vs-model" / "overwash_vs_model.py"
    spec = importlib.util.spec_from_file_location("overwash_vs_model", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# A run's saved model and its metadata
def load_run(run_dir):
    name = run_dir.name
    c = np.load(run_dir / f"{name}.npz", allow_pickle=True)["cascade"][0]
    meta = json.loads((run_dir / f"{name}_run_metadata.json").read_text(encoding="utf-8"))
    return c, meta


# Hit rate against the imagery for one run, by overwash_vs_model's rule
def observed_scores(ovm, c, start, obs):
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    win = (start, 2010 if start == 1996 else 2024)
    q = np.array([np.asarray(b.QowTS) for b in c.barrier3d]).T
    largest, _, _ = ovm.storm_dates(win)
    rows = []
    for _, im in ovm.images_in(win, obs).iterrows():
        o = obs[obs.Obs_ID == im.Obs_ID].set_index("domain").overwash
        for gis in range(1, 91):
            p = DOM.gis_to_pad(gis)
            m = any(q[t, p] > 0 and im["from"] < largest[t] <= im.date
                    for t in range(1, q.shape[0]) if t in largest)
            rows.append(dict(gis=gis, observed=o.get(gis, np.nan), model=int(m)))
    return ovm.scores(pd.DataFrame(rows))


# Fixed against control per member: the tables, then the figure
def compare():
    import HAT_storm_length_selection as S
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    from site_layer import hat_overwash as ow
    ovm = _overwash_module()
    obs = pd.read_csv(ow.OBSERVATIONS / "overwash_observations.csv")
    rp = slice(DOM.gis_to_pad(1), DOM.gis_to_pad(90) + 1)
    rows, per_domain = [], []
    for member, (start, scenario, _) in MEMBERS.items():
        out = {}
        for arm, d in (("control", control_dir(member)), ("fixed", fixed_dir(member))):
            c, meta = load_run(d)
            xs = np.array([np.asarray(b.x_s_TS) for b in c.barrier3d]).T * 10
            q = np.array([np.asarray(b.QowTS) for b in c.barrier3d]).T
            sk = {**S.rescored_skill(d), "lrr_r2_median": meta["skill"]["lrr_r2_median"]}   # today's LOWESS window
            out[arm] = dict(
                net=-(xs[-1] - xs[0])[rp],                      # + seaward, m
                qow=q[1:, rp],
                barrier3d=meta["identity"].get("barrier3d_branch"),
                commit=str(meta["identity"].get("barrier3d_commit"))[:7],
                skill=sk, obs=observed_scores(ovm, c, start, obs),
                years=xs.shape[0] - 1,
            )
        dnet = out["fixed"]["net"] - out["control"]["net"]
        for arm in ("control", "fixed"):
            o = out[arm]
            rows.append(dict(
                member=member, arm=arm, barrier3d=f"{o['barrier3d']}@{o['commit']}",
                mean_net_change_m=float(o["net"].mean()),
                domain_years_overwash=int((o["qow"] > 0).sum()),
                total_overwash_m3_per_m=float(o["qow"].sum()),
                obs_hit_rate=o["obs"]["hit_rate"], obs_both=o["obs"]["both"],
                obs_observed_only=o["obs"]["observed_only"], obs_model_only=o["obs"]["model_only"],
                **{f"skill_{k}": float(v) for k, v in o["skill"].items()
                   if k in ("mean_bias_interior_m_yr", "rmse_interior_m_yr", "lrr_r2_median")},
            ))
        for i, g in enumerate(range(1, 91)):
            per_domain.append(dict(member=member, gis=g, net_control_m=out["control"]["net"][i],
                                   net_fixed_m=out["fixed"]["net"][i], fixed_minus_control_m=dnet[i],
                                   overwash_years_control=int((out["control"]["qow"][:, i] > 0).sum()),
                                   overwash_years_fixed=int((out["fixed"]["qow"][:, i] > 0).sum())))
        print(f"{member:14s} fixed-control net change: mean {dnet.mean():+.2f} m, "
              f"max |.| {np.abs(dnet).max():.1f} m, domains >1 m: {(np.abs(dnet) > 1).sum()}")
    tables = EXP_DIR / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    comp = pd.DataFrame(rows)
    comp.to_csv(tables / "comparison.csv", index=False)
    pd.DataFrame(per_domain).to_csv(tables / "per_domain.csv", index=False)
    with pd.option_context("display.width", 220, "display.max_columns", 40):
        print(comp.round(3).to_string(index=False))
    figure(pd.DataFrame(per_domain))


# Net change and overwash years per domain, control against fixed
def figure(pdm):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer.hat_figure_style import (apply_style, C, INK_MUTED, DOMAIN_AXIS_LABEL,
                                             figsize, save, record_caption, _title, open_frame,
                                             town_bands)
    apply_style()
    fig, axes = plt.subplots(len(MEMBERS), 2, figsize=figsize("double", height=8.0), constrained_layout=True,
                             sharex=True, gridspec_kw=dict(width_ratios=[1.4, 1]))
    for i, member in enumerate(MEMBERS):
        d = pdm[pdm.member == member]
        ax = axes[i, 0]
        ax.plot(d.gis, d.net_control_m, color=C["BASE"], lw=1.2, label="current Barrier3D")
        ax.plot(d.gis, d.net_fixed_m, color=C["ACCENT"], lw=1.2, label="gap + momentum fixes")
        ax.axhline(0, color=INK_MUTED, lw=0.5)
        ax.set_ylim(np.nanpercentile(np.r_[d.net_control_m, d.net_fixed_m], [1, 99]) * 1.3)
        ax.set_ylabel("net change (m)")
        town_bands(ax, label=(i == 0))
        open_frame(ax)
        _title(ax, 2 * i, member.replace("_", " ") + ", shoreline change (+ seaward)")
        ax2 = axes[i, 1]
        ax2.bar(d.gis - 0.2, d.overwash_years_control, width=0.4, color=C["BASE"])
        ax2.bar(d.gis + 0.2, d.overwash_years_fixed, width=0.4, color=C["ACCENT"])
        ax2.set_ylabel("years with overwash")
        open_frame(ax2)
        _title(ax2, 2 * i + 1, "overwash years")
    axes[0, 0].legend(frameon=False, fontsize=7.5, loc="lower left")
    for ax in axes[-1]:
        ax.set_xlabel(DOMAIN_AXIS_LABEL)
    out = save(fig, EXP_DIR / "figures" / "fix_effect.png")
    plt.close(fig)
    record_caption(out[0],
        "What the three Barrier3D overwash fixes (DuneGaps keeps every overtopped cell; gap discharge on "
        "every gap cell; the inundation momentum constant no longer reset to 0) change in the hindcast: "
        "each matrix run (grey, current Barrier3D 49fd069) against its twin on the fix branch (purple, "
        "fix/overwash-gaps-momentum), edgeBE, option A waves. Left: net shoreline change over the window, "
        "positive seaward. Right: the number of model years each domain overwashed. GIS 1 is Cape Point, "
        "GIS 90 Pea Island.")


# Run: the runs or the comparison
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["run", "compare"])
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--only", nargs="*", choices=list(MEMBERS))
    a = ap.parse_args()
    if a.action == "run":
        if not (WORKTREE / "barrier3d" / "barrier3d.py").exists():
            raise SystemExit(f"no Barrier3D worktree at {WORKTREE}")
        with ThreadPoolExecutor(a.workers) as ex:
            list(ex.map(run_member, a.only or list(MEMBERS)))
    else:
        compare()


if __name__ == "__main__":
    main()
