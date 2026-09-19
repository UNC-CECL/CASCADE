"""
HAT_be_dune_edgesolve_results.py
==============================================================================
Close the books on experiments/2026-09-16-dune-edgesolve: which step stands
as each solve's answer, what the pair is, and how the run scores against BOTH
observations.

WHAT IT WRITES  (under output/raw_runs/experiments/2026-09-16-dune-edgesolve/)
    solved.csv      one row per (window, smooth): the solved step, its run
                    name, the pair at GIS 1 / 90, and the CoastSat-solved pair
                    it replaces. HAT_rate_windows.py reads this to find the
                    dune-solved runs (model sets dune-mean3, dune-raw).
    skill.csv       every solve AND its CoastSat-solved counterpart scored the
                    same two ways over GIS 2-89: against the CoastSat LRR
                    target (model lrr_m_yr, as run_index.csv scores) and
                    against the dune-line endpoint rate (model
                    change_rate_m_yr, as HAT_rate_windows.py vs_duneline/net_change
                    scores).
    RESULTS.md      the two tables, rendered.

USAGE
    python HAT_be_dune_edgesolve_results.py --solved 1984:raw:3 1984:mean3:2 ...
        # window start year : smoothing : the step that converged
    python HAT_be_dune_edgesolve_results.py --exp 2026-09-18-dune-edgesolve         --solved 1984:raw:3@2026-09-16-dune-edgesolve 1996:raw:2 ...
        # --exp is where the files are written and where a bare spec's runs
        # sit; "@<experiment>" carries a solve over from an earlier one (the
        # 09-18 re-solve kept 1984-2004, whose lines did not change)

    The dune target is read from 5-scr/3-rates/duneline/endpoint/ (2026-09-18),
    the same stored product HAT_rate_windows.py draws.

Author: Hannah A. Henry, UNC CECL
==============================================================================
"""
from __future__ import annotations

import argparse
import datetime as _dt
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
PROJECT_ROOT = next(_p for _p in _HERE.parents if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
sys.path.insert(0, str(PROJECT_ROOT))

from site_layer.hatteras_site_config import HATTERAS_PERIODS, HATTERAS_BE_EDGE_DOMAINS  # noqa: E402
from cascade_pipeline.run_registry import find_run_dir, load_run_index    # noqa: E402

RUN_ROOT = PROJECT_ROOT / "output" / "raw_runs"
EXP = "2026-09-16-dune-edgesolve"          # default; --exp overrides
EXP_DIR = RUN_ROOT / "experiments" / EXP
# The 1984 and 2004 brackets were run once, under the 09-16 experiment; their
# inputs did not change with the 09-18 re-digitization (the 1984 and 2004
# lines, and the run itself never reads a dune line), so every re-solve
# reuses them from there.
BRACKET_EXP = "2026-09-16-dune-edgesolve"
INTERIOR = (2, 89)

# The CoastSat-solved counterpart of each window: the matrix run for 1996
# and 2010, the fresh bracket for 1984 and 2004 (same code, same versions).
RUN_NAME = {
    1984: "HAT_1984_2004_edgeBE_road_bdm_nogroin",
    1996: "HAT_1996_2010_edgeBE_road_bdm_nogroin",
    2004: "HAT_2004_2024_edgeBE_road_bdm_nourish_nogroin",
    2010: "HAT_2010_2024_edgeBE_road_bdm_nourish_nogroin",
}
COASTSAT_RUN = {
    1984: ("experiment", f"{BRACKET_EXP}/brackets"),
    1996: ("matrix", ""),
    2004: ("experiment", f"{BRACKET_EXP}/brackets"),
    2010: ("matrix", ""),
}


def _solve_module():
    path = _HERE.with_name("HAT_be_edge_domain_solve.py")
    spec = importlib.util.spec_from_file_location("HAT_be_edge_domain_solve", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def rates_table(start, kind, tag):
    end = HATTERAS_PERIODS[start]["end_year"]
    run_dir = find_run_dir(RUN_ROOT, RUN_NAME[start], (start, end), "edgeBE",
                           kind=kind, tag=tag)
    return pd.read_csv(run_dir / "tables" / "shoreline_change_rate.csv").set_index("gis_domain")


def index_row(index, start, kind, tag):
    rows = index[(index["run_name"] == RUN_NAME[start]) & (index["kind"] == kind)
                 & (index["tag"] == tag)]
    if rows.empty:
        raise ValueError(f"{RUN_NAME[start]} ({kind}, {tag!r}) not in run_index.csv")
    row = rows.iloc[-1].copy()
    # the index loader keeps every column as text; the two rates are numbers here
    for col in ("be_rate_gis1_m_yr", "be_rate_gis90_m_yr"):
        row[col] = float(row[col]) if str(row[col]).strip() not in ("", "nan") else 0.0
    return row


def score(model, target, column):
    lo, hi = INTERIOR
    ids = [g for g in range(lo, hi + 1) if g in target and not np.isnan(target[g])]
    r = np.array([model.loc[g, column] - target[g] for g in ids])
    return float(r.mean()), float(np.sqrt((r ** 2).mean())), len(ids)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--solved", nargs="+", required=True,
                    help="start:smooth:step[@experiment], e.g. 1984:raw:3")
    ap.add_argument("--exp", default=EXP,
                    help="the experiment the files are written to")
    args = ap.parse_args(argv)
    global EXP_DIR
    EXP_DIR = RUN_ROOT / "experiments" / args.exp
    solve = _solve_module()
    index = load_run_index(RUN_ROOT / "run_index.csv")

    solved_rows, skill_rows = [], []
    targets = {}
    for start in sorted(RUN_NAME):
        end = HATTERAS_PERIODS[start]["end_year"]
        cs_target = solve.load_target(start, end)
        # full per-domain dune rate, raw (the interior score does not smooth),
        # from the stored product
        from site_layer.hat_observed_rates import dune_endpoint_csv
        rate = (pd.read_csv(dune_endpoint_csv(start, end, "domain"))
                .set_index("domain_number")["mean_rate_m_yr"])
        dune_all = {int(g): float(v) for g, v in rate.items()}
        targets[start] = (cs_target, dune_all)

        # the CoastSat-solved counterpart
        kind, tag = COASTSAT_RUN[start]
        row = index_row(index, start, kind, tag)
        model = rates_table(start, kind, tag)
        for obs, target, col in (("coastsat_lrr", cs_target, "lrr_m_yr"),
                                 ("duneline_endpoint", dune_all, "change_rate_m_yr")):
            b, e, n = score(model, target, col)
            skill_rows.append(dict(window=f"{start}_{end}", solve="coastsat", step="",
                                   gis1=row["be_rate_gis1_m_yr"], gis90=row["be_rate_gis90_m_yr"],
                                   scored_against=obs, model_column=col, n=n,
                                   bias_m_yr=b, rmse_m_yr=e))

    for spec in args.solved:
        spec, _, from_exp = spec.partition("@")
        start_s, smooth, step_s = spec.split(":")
        start, step = int(start_s), int(step_s)
        end = HATTERAS_PERIODS[start]["end_year"]
        tag = f"{from_exp or args.exp}/{smooth}/step{step}"
        row = index_row(index, start, "experiment", tag)
        model = rates_table(start, "experiment", tag)
        cs_kind, cs_tag = COASTSAT_RUN[start]
        cs_row = index_row(index, start, cs_kind, cs_tag)
        dune_target, _ = solve.load_dune_target(start, end, smooth)
        solved_rows.append(dict(
            window=f"{start}_{end}", smooth=smooth, step=step, run_name=RUN_NAME[start],
            tag=tag, gis1=row["be_rate_gis1_m_yr"], gis90=row["be_rate_gis90_m_yr"],
            target_gis1=dune_target[HATTERAS_BE_EDGE_DOMAINS[0]],
            target_gis90=dune_target[HATTERAS_BE_EDGE_DOMAINS[1]],
            model_gis1=model.loc[HATTERAS_BE_EDGE_DOMAINS[0], "change_rate_m_yr"],
            model_gis90=model.loc[HATTERAS_BE_EDGE_DOMAINS[1], "change_rate_m_yr"],
            coastsat_gis1=cs_row["be_rate_gis1_m_yr"], coastsat_gis90=cs_row["be_rate_gis90_m_yr"],
            timestamp=row["timestamp"], git_commit=row["git_commit"]))
        cs_target, dune_all = targets[start]
        for obs, target, col in (("coastsat_lrr", cs_target, "lrr_m_yr"),
                                 ("duneline_endpoint", dune_all, "change_rate_m_yr")):
            b, e, n = score(model, target, col)
            skill_rows.append(dict(window=f"{start}_{end}", solve=f"dune-{smooth}", step=step,
                                   gis1=row["be_rate_gis1_m_yr"], gis90=row["be_rate_gis90_m_yr"],
                                   scored_against=obs, model_column=col, n=n,
                                   bias_m_yr=b, rmse_m_yr=e))

    solved = pd.DataFrame(solved_rows).sort_values(["window", "smooth"])
    skill = pd.DataFrame(skill_rows).sort_values(["window", "solve", "scored_against"])
    solved.to_csv(EXP_DIR / "solved.csv", index=False)
    skill.to_csv(EXP_DIR / "skill.csv", index=False)

    # RESULTS.md
    f = lambda v: f"{v:+.1f}"  # noqa: E731
    g = lambda v: f"{v:+.2f}"  # noqa: E731
    lines = [f"# {args.exp} - results", "",
             f"Written {_dt.date.today().isoformat()} by "
             "scripts/input_prep/7-source-sink/2-calibrate/HAT_be_dune_edgesolve_results.py. "
             "See NOTE.md for the question and the layout.", "",
             "## The solved pairs, GIS 1 / GIS 90, m/yr", "",
             "| window | reading | step | dune-line solve | target at the ends | CoastSat solve |",
             "|---|---|---|---|---|---|"]
    for _, r in solved.iterrows():
        lines.append(f"| {r.window.replace('_', '-')} | {r.smooth} | {r.step} | "
                     f"{f(r.gis1)} / {f(r.gis90)} | {g(r.target_gis1)} / {g(r.target_gis90)} | "
                     f"{f(r.coastsat_gis1)} / {f(r.coastsat_gis90)} |")
    lines += ["", "## Interior skill, GIS 2-89, model minus observation, m/yr", "",
              "Each run scored against both observations: the CoastSat LRR target "
              "(model OLS slope, as run_index.csv) and the dune-line endpoint rate "
              "(model endpoint rate, as model_vs_observed/vs_duneline/net_change). The CoastSat "
              "row is the run the dune solve started from.", "",
              "| window | solve | ends | vs CoastSat bias | RMSE | vs dune line bias | RMSE |",
              "|---|---|---|---|---|---|---|"]
    for (w, sv), grp in skill.groupby(["window", "solve"], sort=False):
        cs = grp[grp.scored_against == "coastsat_lrr"].iloc[0]
        dl = grp[grp.scored_against == "duneline_endpoint"].iloc[0]
        lines.append(f"| {w.replace('_', '-')} | {sv} | {f(cs.gis1)} / {f(cs.gis90)} | "
                     f"{g(cs.bias_m_yr)} | {cs.rmse_m_yr:.2f} | {g(dl.bias_m_yr)} | {dl.rmse_m_yr:.2f} |")
    (EXP_DIR / "RESULTS.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
