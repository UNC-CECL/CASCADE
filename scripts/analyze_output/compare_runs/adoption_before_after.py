r"""
adoption_before_after.py -- the matrix before and after the 2026-09-28 adoption
==============================================================================
Hannah, 2026-09-28: "show me the comparison when the matrix finishes".

    before  output/raw_runs/archive/2026-09-28-pre-ceiling/matrix/: Barrier3D
            fix/route-overwash-axis-swap (49fd069), Dmaxel default (3.4 m
            NAVD88), storms v3_72, the pre-adoption LOWESS-7 ends
    after   output/raw_runs/matrix/: Barrier3D hatteras/adopted (overwash
            fixes + per-cell dune ceilings), storms v3_trim24, the ends
            re-solved on it (end-domain-boundaries/2026-09-28-ends-resolved-adopted)

For every matrix run (both windows, both presets, every scenario):
    shoreline  interior (GIS 2-89) RMSE and bias of the model LRR against the
               CoastSat LOWESS-7 target (run_registry.skill_vs_target, as the
               runner scores), and the spatial correlation r
    overwash   against the imagery (8-overwash-analysis), each run dated by its
               own storm file: POD, POFD, PSS, timing r, space r
    dunes      1996-2010 runs: the 2010 dune crest minus the 2009 lidar

Per-domain tables beside the scores: cells_<side>.csv (every image x domain,
observed and model overwash) and crest_<side>.csv (end-of-run crest per GIS
domain, m MHW); crest_lidar_2009.csv is the 2010-start dune file's crest.
Figures: adoption_before_after_figures.py.

Each side is scored under the Barrier3D it ran on, so the storm sharing uses
that version's DuneGaps and DuneGrowth: `score --side before` must run with
PYTHONPATH=<Barrier3D at 49fd069 + the ceiling feature, off> (the worktree
../Barrier3D-dune-ceiling), `score --side after` with the editable install.

    python adoption_before_after.py score --side before   (PYTHONPATH=../Barrier3D-dune-ceiling)
    python adoption_before_after.py score --side after
    python adoption_before_after.py report

WHERE: output/comparisons/adoption_2026-09-28/
==============================================================================

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
PROJECT_ROOT = next(_p for _p in _HERE.parents if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "hatteras_ms" / "experiments"))
import HAT_storm_length_selection as S  # noqa: E402

RAW = PROJECT_ROOT / "output" / "raw_runs"
SIDES = {"before": RAW / "archive" / "2026-09-28-pre-ceiling" / "matrix", "after": RAW / "matrix"}
STORM_VARIANT = {"before": "v3_72", "after": "v3_trim24"}
OUT = PROJECT_ROOT / "output" / "comparisons" / "adoption_2026-09-28"
WINDOWS = ((1996, 2010), (2010, 2024))


def runs(side):
    for w in WINDOWS:
        for preset_dir in sorted((SIDES[side] / S.wtag(w)).glob("*")):
            for rd in sorted(preset_dir.glob("HAT_*")):
                if (rd / f"{rd.name}.npz").exists():
                    yield w, preset_dir.name, rd


def target_table(w, _cache={}):
    if w not in _cache:
        from cascade_pipeline.hindcast import build_target_table
        from cascade_pipeline.coastsat_lowess import CoastSatDataset, LowessConfig, build_coastsat_series
        from site_layer.hatteras_site_config import HATTERAS_DOMAINS
        from site_layer.hat_observed_rates import lrr_csv
        cfg = LowessConfig(window_domains=(7,), skip_southern_domains=10)
        series = build_coastsat_series([CoastSatDataset(label=f"CoastSat {w[0]}", period_start=w[0],
                                                        csv_path=str(lrr_csv(*w)))], w[0], cfg,
                                       domains=HATTERAS_DOMAINS)
        _cache[w] = build_target_table(series[0], cfg, HATTERAS_DOMAINS, 7)
    return _cache[w]


def shoreline(rd, w):
    from cascade_pipeline.run_registry import skill_vs_target
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM, SCORE_INTERIOR_GIS
    lrr = pd.read_csv(rd / "tables" / "shoreline_change_rate.csv").set_index("gis_domain")["lrr_m_yr"]
    model = np.full(DOM.total_domains, np.nan)
    for g, v in lrr.items():
        model[DOM.gis_to_pad(int(g))] = v
    tab = target_table(w)
    sk = skill_vs_target(model, tab, DOM, interior_gis=SCORE_INTERIOR_GIS)
    lo, hi = SCORE_INTERIOR_GIS
    t = tab.set_index("gis_domain")["target_lrr_m_yr"]
    g = [x for x in range(lo, hi + 1) if x in t.index and x in lrr.index]
    r = float(np.corrcoef(lrr[g].values, t[g].values)[0, 1])
    return dict(rmse_interior_m_yr=float(sk["rmse_interior_m_yr"]),
                bias_interior_m_yr=float(sk["mean_bias_interior_m_yr"]), shoreline_r=r)


def score(side):
    import barrier3d
    import barrier3d.barrier3d as b3d
    import inspect
    src = inspect.getsource(b3d)
    has_fixes = "start:stop + 1" in src
    if side == "before" and has_fixes:
        raise SystemExit("score the BEFORE side under the pre-adoption DuneGaps: "
                         "PYTHONPATH=../Barrier3D-dune-ceiling")
    if side == "after" and not has_fixes:
        raise SystemExit("score the AFTER side under hatteras/adopted (the editable install)")
    from site_layer import hat_overwash as ow, hat_env_forcings as env
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    ovm = S.overwash_module()
    obs = pd.read_csv(ow.OBSERVATIONS / "overwash_observations.csv")
    pads = [DOM.gis_to_pad(g) for g in range(1, 91)]
    lid = None
    rows, all_cells, crests = [], [], []
    for w, preset, rd in runs(side):
        c = S.load_state(rd)
        meta = json.loads(next(rd.glob("*_run_metadata.json")).read_text(encoding="utf-8"))
        summ = pd.read_csv(env.storm_summary_file(*w, variant=STORM_VARIANT[side]))
        cells = S.overwash_cells(c, w, summ, obs, ovm)
        row = dict(side=side, window=S.wtag(w), preset=preset, run=rd.name,
                   scenario=meta["index row"].get("scenario"),
                   barrier3d=f"{meta['identity'].get('barrier3d_branch')}@{str(meta['identity'].get('barrier3d_commit'))[:7]}",
                   be_gis1=meta["index row"].get("be_rate_gis1_m_yr"), be_gis90=meta["index row"].get("be_rate_gis90_m_yr"),
                   **shoreline(rd, w), **S.overwash_scores(cells, 0.0))
        crest = np.array([np.median((np.asarray(c.barrier3d[p].DuneDomain[-1]).max(axis=1)
                                     + c.barrier3d[p].BermEl) * 10) for p in pads])
        row["crest_end_median_m"] = float(np.median(crest))
        all_cells.append(cells.assign(window=S.wtag(w), preset=preset, run=rd.name))
        crests.append(pd.DataFrame(dict(window=S.wtag(w), preset=preset, run=rd.name,
                                        gis=range(1, 91), crest_end_m_mhw=crest)))
        if w[0] == 1996:
            if lid is None:
                lc = S.load_state(next(r for ww, pp, r in runs("after") if ww[0] == 2010))
                lid = np.array([np.median((np.asarray(lc.barrier3d[p].DuneDomain[0]).max(axis=1)
                                           + lc.barrier3d[p].BermEl) * 10) for p in pads])
            row["crest_2010_minus_lidar_m"] = float(np.median(crest - lid))
        rows.append(row)
        print(f"  {side} {S.wtag(w)} {preset} {rd.name}", flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUT / f"scores_{side}.csv", index=False)
    pd.concat(all_cells).to_csv(OUT / f"cells_{side}.csv", index=False)
    pd.concat(crests).to_csv(OUT / f"crest_{side}.csv", index=False)
    if side == "after":
        pd.DataFrame(dict(gis=range(1, 91), crest_2009_lidar_m_mhw=lid)).to_csv(OUT / "crest_lidar_2009.csv", index=False)


def report():
    b = pd.read_csv(OUT / "scores_before.csv")
    a = pd.read_csv(OUT / "scores_after.csv")
    key = ["window", "preset", "run"]
    m = b.merge(a, on=key, suffixes=("_before", "_after"), how="outer")
    cols = ["rmse_interior_m_yr", "bias_interior_m_yr", "shoreline_r", "PSS", "POD", "POFD", "timing_r", "space_r"]
    out = m[key + [f"{c}_{s}" for c in cols for s in ("before", "after")]]
    out.to_csv(OUT / "before_after.csv", index=False)
    short = m.copy()
    short["scenario"] = short["run"].str.replace(r"HAT_\d{4}_\d{4}_\w+?_offsetmetres_", "", regex=True)
    with pd.option_context("display.width", 260, "display.max_columns", 30):
        for c in cols:
            short[c] = short[f"{c}_before"].round(2).astype(str) + " → " + short[f"{c}_after"].round(2).astype(str)
        print(short[["window", "preset", "scenario"] + cols].to_string(index=False))
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["score", "report"])
    ap.add_argument("--side", choices=list(SIDES))
    a = ap.parse_args()
    score(a.side) if a.action == "score" else report()


if __name__ == "__main__":
    main()
