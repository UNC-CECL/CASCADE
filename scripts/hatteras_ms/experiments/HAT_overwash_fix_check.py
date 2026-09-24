r"""
HAT_overwash_fix_check.py -- how much does the Barrier3D route_overwash fix move the results?
==============================================================================
THE BUG (found 2026-09-24): Barrier3D/barrier3d/barrier3d.py, route_overwash,
the subaerial test indexed Elevation[TS, i, d+1:d+10] (row i, columns d+1..d+9)
where Elevation[TS, d+1:d+10, i] (the nine cells landward of the flow) is
meant: the wrong cells whenever i < rows, out of bounds whenever i >= rows.
The fix is one line, commit 49fd069 on Barrier3D branch
fix/route-overwash-axis-swap (local). Barrier3D is installed editable, so the
branch checked out IS the model every run uses.

THE CHECK (Hannah: "patch it on a branch and measure the impact")
    patched            8 runs on the fix branch, each the twin of a run
                       already made unpatched
    patched_boundscheck the natural 1996 and managed 2010 baselines again with
                       NUMBA_BOUNDSCHECK=1: does the fixed model read out of
                       bounds anywhere else?
    unpatched          the two /10 twins re-run on master with today's code:
                       their archived matrix runs were made on older CASCADE
                       code, so they are not a clean control. The six metres
                       twins were made today with today's code and are.
    compare            runner scores, share of variation explained, and the
                       per-domain LRR difference, patched minus unpatched

    Every launch records the Barrier3D branch and commit it ran on, and
    refuses to run a patched member off the fix branch or an unpatched one on it.

WHERE: output/raw_runs/experiments/2026-09-24-overwash-fix/
           <member>/<period>/<preset>/<run_name>/   runs (on disk only)
           logs/, launches.jsonl, comparison.csv, NOTE.md

USAGE
    python HAT_overwash_fix_check.py run patched            (fix branch checked out)
    python HAT_overwash_fix_check.py run unpatched          (master checked out)
    python HAT_overwash_fix_check.py compare
==============================================================================
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve()
PROJECT_ROOT = next(_p for _p in _HERE.parents if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "sensitivity_analysis"))
sys.path.insert(0, str(_HERE.parent))
import natural_wave_sensitivity as N  # noqa: E402
import HAT_offset_scale_wave_tuning as common  # noqa: E402

BARRIER3D = PROJECT_ROOT.parent / "Barrier3D"
FIX_BRANCH = "fix/route-overwash-axis-swap"
TAG = "2026-09-24-overwash-fix"
EXP_DIR = PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
RAW = PROJECT_ROOT / "output" / "raw_runs"
STUDY = N.STUDY_DIR
ARCH = RAW / "archive" / "2026-09-24-pre-metres" / "matrix"
DIV10 = {"HAT_OFFSET_MODE": "asrun"}          # plus the build the /10 runs used, per period

# member -> (period, scenario, wave settings or None for the model defaults,
#            extra env, the unpatched twin's run folder or None)
BASE = dict(N.BASELINE)
MEMBERS = {
    "natural_baseline_1996": (1996, N.SCENARIO, BASE, {},
                              STUDY / "baseline/1996_2010/zeroBE"),
    "managed_baseline_1996": (1996, N.MANAGED, BASE, {},
                              STUDY / "baseline_full_management/1996_2010/zeroBE"),
    "natural_baseline_2010": (2010, N.SCENARIO, BASE, {},
                              STUDY / "baseline/2010_2024/zeroBE"),
    "managed_baseline_2010": (2010, N.MANAGED, BASE, {},
                              STUDY / "baseline_full_management/2010_2024/zeroBE"),
    "natural_ridge_1996":    (1996, N.SCENARIO,
                              {**BASE, "hs": 1.25, "wave_angle_high_fraction": 0.5}, {},
                              STUDY / "grid_wave_height_x_high_angle/1996_2010/zeroBE"),
    "managed_highangle0.4_2010": (2010, N.MANAGED, {**BASE, "wave_angle_high_fraction": 0.4},
                                  {}, None),   # crashed unpatched in year 13
    "div10_managed_1996":    (1996, N.MANAGED, None,
                              {**DIV10, "HAT_OFFSET_VERSION_1996": "superseded_20260924_pre-metres/v1"},
                              "unpatched"),
    "div10_managed_2010":    (2010, N.MANAGED, None,
                              {**DIV10, "HAT_OFFSET_VERSION_2010": "superseded_20260924_pre-metres/v1"},
                              "unpatched"),
}
BOUNDSCHECK = ("natural_baseline_1996", "managed_baseline_2010")


def barrier3d_state():
    def git(*a):
        return subprocess.run(["git", "-C", str(BARRIER3D), *a], capture_output=True,
                              text=True).stdout.strip()
    return {"branch": git("branch", "--show-current"), "commit": git("rev-parse", "--short", "HEAD"),
            "dirty": bool(git("status", "--short", "barrier3d"))}


def env_for(member, variant):
    period, scenario, waves, extra, _ = MEMBERS[member]
    env = N.run_env("x", period, dict(BASE), scenario)
    if waves is None:                             # /10 twins: the model's own wave defaults
        for k in N.ENV.values():
            env.pop(k, None)
        env["HAT_OFFSET_MODE"] = "asrun"
    else:
        env.update({N.ENV[k]: f"{v}" for k, v in waves.items()})
    env.update(extra)
    env["HAT_RUN_KIND"] = "experiment"
    env["HAT_RUN_TAG"] = f"{TAG}/{variant}_{member}"
    if variant == "patched_boundscheck":
        env["NUMBA_BOUNDSCHECK"] = "1"
    return env


def launch(member, variant):
    state = barrier3d_state()
    want_fix = variant.startswith("patched")
    if (state["branch"] == FIX_BRANCH) != want_fix or state["dirty"]:
        raise SystemExit(f"{variant} needs Barrier3D on "
                         f"{FIX_BRANCH if want_fix else 'master'} and clean; it is on "
                         f"{state['branch']} (dirty={state['dirty']})")
    log = EXP_DIR / "logs" / f"{variant}_{member}.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    p = subprocess.run([sys.executable, str(N.HINDCAST)], env=env_for(member, variant),
                       cwd=str(PROJECT_ROOT), capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=3600)
    log.write_text((p.stdout or "") + "\n--- STDERR ---\n" + (p.stderr or ""), encoding="utf-8")
    rec = dict(member=member, variant=variant, barrier3d=state, exit=p.returncode,
               minutes=round((time.perf_counter() - t0) / 60, 1),
               outcome="finished" if p.returncode == 0 else N.stop_reason(log),
               index_error="IndexError" in (p.stderr or ""))
    with (EXP_DIR / "launches.jsonl").open("a", encoding="utf-8") as f:
        f.write(json.dumps(rec) + "\n")
    print(f"{variant} {member}: {rec['outcome']} ({rec['minutes']} min, Barrier3D "
          f"{state['branch']}@{state['commit']})", flush=True)
    return rec


def cmd_run(a):
    common.keep_awake()
    if a.what == "patched":
        jobs = [(m, "patched") for m in MEMBERS] + [(m, "patched_boundscheck") for m in BOUNDSCHECK]
    else:
        jobs = [(m, "unpatched") for m, v in MEMBERS.items() if v[4] == "unpatched"]
    with ThreadPoolExecutor(max_workers=a.jobs) as pool:
        list(pool.map(lambda j: launch(*j), jobs))
    return 0


def run_dir(member, variant):
    """The one run folder under experiments/<TAG>/<variant>_<member>/."""
    hits = sorted((EXP_DIR / f"{variant}_{member}").glob("*/*/*/*_run_metadata.json"))
    return hits[0].parent if hits else None


def twin(member):
    period, scenario, waves, extra, ref = MEMBERS[member]
    if ref is None:
        return None
    if ref == "unpatched":
        return run_dir(member, "unpatched")
    for d in sorted(ref.iterdir()):
        md = json.loads(next(d.glob("*_run_metadata.json")).read_text(encoding="utf-8")) \
            if list(d.glob("*_run_metadata.json")) else None
        if md is None:
            continue
        w = md["wave climate"]
        got = {"hs": float(w["wave_height_m"]), "wave_period_s": float(w["wave_period_s"]),
               "wave_asymmetry": float(w["wave_asymmetry"]),
               "wave_angle_high_fraction": float(w["wave_angle_high_frac"])}
        if all(np.isclose(got[k], waves[k]) for k in waves):
            return d
    return None


def scores(d, period):
    rates = common.run_rates(d)
    sc = common.alongshore_scores(rates, common.coastsat_target(period))
    idx = json.loads(next(d.glob("*_run_metadata.json")).read_text(encoding="utf-8"))["index row"]
    return rates, dict(bias=float(idx["mean_bias_interior_m_yr"]),
                       rmse=float(idx["rmse_interior_m_yr"]), ve=sc["variance_explained"],
                       r=sc["r_alongshore"])


def cmd_compare(a):
    import pandas as pd
    launches = [json.loads(l) for l in (EXP_DIR / "launches.jsonl").read_text().splitlines()]
    rows = []
    for member, (period, scenario, *_rest) in MEMBERS.items():
        row = dict(member=member, period=f"{period}-{period + 14}", scenario=scenario)
        pd_ = run_dir(member, "patched")
        tw = twin(member)
        row["patched"] = "scored" if pd_ else next((l["outcome"] for l in launches
                                                    if l["member"] == member and l["variant"] == "patched"), "not run")
        row["unpatched"] = "scored" if tw else "crashed (year 13)" if MEMBERS[member][4] is None else "missing"
        if pd_ is not None:
            rp, sp = scores(pd_, period)
            row.update({f"patched_{k}": v for k, v in sp.items()})
        if tw is not None:
            ru, su = scores(tw, period)
            row.update({f"unpatched_{k}": v for k, v in su.items()})
        if pd_ is not None and tw is not None:
            d = common.interior(rp) - common.interior(ru)
            row.update(rate_change_rms=float(np.sqrt((d ** 2).mean())),
                       rate_change_max_abs=float(d.abs().max()),
                       rate_change_mean=float(d.mean()),
                       domains_moved_over_0p1=int((d.abs() > 0.1).sum()),
                       domains_moved_over_0p5=int((d.abs() > 0.5).sum()),
                       run_patched=str(pd_.relative_to(RAW)), run_unpatched=str(tw.relative_to(RAW)))
        rows.append(row)
    for l in launches:
        if l["variant"] == "patched_boundscheck":
            rows.append(dict(member=l["member"], period="", scenario="bounds check on the fix",
                             patched=l["outcome"],
                             unpatched="IndexError (earlier today)"))
    out = pd.DataFrame(rows)
    out.to_csv(EXP_DIR / "comparison.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 30,
                           "display.float_format", "{:.3f}".format):
        print(out.drop(columns=[c for c in out.columns if c.startswith("run_")]).to_string(index=False))
    return 0


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("what", choices=("patched", "unpatched"))
    r.add_argument("--jobs", type=int, default=6)
    sub.add_parser("compare")
    a = p.parse_args()
    return cmd_run(a) if a.cmd == "run" else cmd_compare(a)


if __name__ == "__main__":
    sys.exit(main())
