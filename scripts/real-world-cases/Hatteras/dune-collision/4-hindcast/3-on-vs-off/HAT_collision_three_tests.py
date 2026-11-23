"""
A simple scorecard for the collision designs (2026-11-19): three yes/no tests on every design and Cs of the sweep.

    python scripts/real-world-cases/Hatteras/dune-collision/4-hindcast/3-on-vs-off/HAT_collision_three_tests.py lidar
    python scripts/real-world-cases/Hatteras/dune-collision/4-hindcast/3-on-vs-off/HAT_collision_three_tests.py score

1. Beats no collision: island crest RMSE (HAT_collision_design_cs_sweep.score) below no collision's on the same
   starting rows.
2. Stable: the crest RMSE moves less than STABLE_M to the next Cs of the sweep.
3. Dune volume: the dune's sand volume change since 1996 (m3 per m of shoreline) against the lidar, its RMSE below
   no collision's. Added because the crest score reads only the higher dune row (dune_crest_comparison.model_crests:
   max over rows), so a design that flattens the front row and leaves the back row alone scores like no collision
   (the back-row check, 2026-11-19). Volume sees both rows and is what the module removes.

Volumes. Model: the two dune rows' heights above the berm (never below 0), times the 10 m row width, per 10 m cell.
Lidar: per 10 m cell, the ground above the berm (1.7 m NAVD88) summed over the survey's own dune window (the
start's window moved as in survey_dune_crests.csv, clipped to the survey's beach start, as the crests are), times
10 m; only where the whole window was surveyed and the crest was not a spike. The two volumes differ in width (the
model's dune is 20 m, a picked window is often wider), so each is compared as its change since 1996: the model
from its starting state, the lidar from the 1996 survey (the same coarse grid as 1997-2005, so the grids' smoothing
largely cancels; volume is less affected by smoothing than a crest). Per GIS domain and survey: the median over
the cells both have. Scored over GIS 6-90, surveys after 1996, without Buxton's rebuilt-dune cells (GIS 9-14 in
2004, 2005, 2009), as the crest scores are.

"lidar" writes data/hatteras_init/1-barrier3d-domains/1996-start/survey-crests/survey_dune_volumes.csv;
"score" writes 3-collision-design/0-all-designs/supporting/three_tests_1996_2009.csv and prints the scorecard.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-19
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
_HAT = next(p for p in _HERE.parents if (p / "hatteras_paths.py").exists())
sys.path.insert(0, str(_HAT))
sys.path.insert(0, str(_HERE.parent))
sys.path.insert(0, str(_HAT / "input_prep" / "1-barrier3d-domains" / "2-survey-crests"))

# --- CONFIG ------------------------------------------------------------------
START = "1996-start"
GIS_LO, GIS_HI = 6, 90
STABLE_M = 0.10                 # test 2: largest crest-RMSE change to the next Cs
ROW_WIDTH_M = 10.0              # Barrier3D's dune cell (two rows, 20 m)
CELL_M = 10.0                   # lidar profile cell after the 10 m resample
VOLUMES = "survey_dune_volumes.csv"
OUT = "three_tests_1996_2009.csv"
# -----------------------------------------------------------------------------


def lidar():
    """Dune volume above the berm per survey, domain and cell, in the survey windows the crests use."""
    import rasterio
    import HAT_survey_dune_crests as S
    X, tv = S.X, S.tv
    prod, surveys = S.STARTS[START]
    windows = json.loads(tv.picks_dir(START).joinpath(
        f"HAT_dune_search_windows_{tv.topo_dirs(START)[2]}.json").read_text(encoding="utf-8"))
    frames = {g: S.frame(START, g, windows) for g in range(1, 91)}
    grids = {g: S.clip_grid(S.product(prod).gapfill_1m / f"clip_domain_{g}_filled.tif") for g in range(1, 91)}
    crest = pd.read_csv(tv.product_dir(START) / "survey-crests" / "survey_dune_crests.csv",
                        usecols=["survey", "gis", "cell", "status"])
    spike = set(map(tuple, crest[crest.status == "spike"][["survey", "gis", "cell"]].to_numpy()))
    berm = X.BERM_ELEV_NAVD_M - X.MHW_M                                  # m MHW
    rows = []
    for year, (kind, path, *_rest) in surveys.items():
        files = sorted(path.parent.glob(path.name)) if "*" in path.name else [path]
        srcs = [rasterio.open(f) for f in files]
        try:
            unit = "m" if kind == "mosaic" else S.unit_of(srcs, START, grids)[0]
            scale, shift = (S.US_FT if unit == "ft" else 1.0), S.SURVEY_VERTICAL_OFFSET_M.get(START, {}).get(year, 0.0)
            for g in range(1, 91):
                shear, c0, i0, i1 = frames[g]
                on = S.mosaic_on_clip(srcs[0], grids[g]) if kind == "mosaic" else S.raster_on_clip(srcs, grids[g])
                z = S.profiles(S.downsample(on * scale + shift), shear, c0)
                i0s, i1s = S.survey_windows(START, year, g, i0, i1, n=z.shape[0])
                above = z > X.BEACH_START_THR_M
                sb = np.where(above.any(axis=1), above.argmax(axis=1), -1)
                for c in range(min(50, z.shape[0])):
                    a = max(int(i0s[c]), 0)
                    a = max(a, int(sb[c])) if (X.CLIP_WINDOW_TO_BEACH and sb[c] >= 0) else a
                    b = min(int(i1s[c]), z.shape[1])
                    win = z[c, a:b]
                    ok = b > a and not (win <= X.NODATA_SENTINEL_M + 1e-9).any() and (year, g, c) not in spike
                    hgt = np.clip(np.where(win > X.SENTINEL_WATER_M + 1e-9, win, berm) - berm, 0, None)
                    v = float(hgt.sum() * CELL_M) if ok else np.nan
                    # seaward and landward halves of the window (the middle cell left out when it is odd):
                    # mean height above the berm, m, to set against the model's front and back rows
                    k = hgt.size // 2
                    sea = float(hgt[:k].mean()) if ok and k else np.nan
                    land = float(hgt[-k:].mean()) if ok and k else np.nan
                    rows.append(dict(survey=year, gis=g, cell=c, volume_m3_per_m=v, window_cells=max(b - a, 0),
                                     seaward_half_m=sea, landward_half_m=land))
        finally:
            for s in srcs:
                s.close()
        print(f"{year}: volume in {sum(np.isfinite(r['volume_m3_per_m']) for r in rows if r['survey'] == year)} cells",
              flush=True)
    out = tv.product_dir(START) / "survey-crests" / VOLUMES
    pd.DataFrame(rows).round(2).to_csv(out, index=False)
    print(out)


def model_volume_change(v, survey_years):
    """Per survey year: (90 GIS x 50 cells) model dune volume change since its start, m3/m."""
    from cascade_pipeline.saved_runs import load_saved_cascade
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    import HAT_collision_design_cs_sweep as DS
    rd = DS.H.run_dir(v)
    meta = json.loads(next(rd.glob("*_run_metadata.json")).read_text(encoding="utf-8"))
    first = int(meta["period"]["first_model_year"])
    c = load_saved_cascade(rd / f"{rd.name}.npz")
    d = DOM.for_total(len(c.barrier3d))
    vol = np.stack([np.clip(np.asarray(c.barrier3d[d.gis_to_pad(g)].DuneDomain, float) * 10, 0, None).sum(axis=2)
                    * ROW_WIDTH_M for g in range(1, 91)], axis=1)                       # states x 90 x 50
    return {y: vol[s - first] - vol[0] for y, s in survey_years.items() if 0 <= s - first < vol.shape[0]}


def score():
    import HAT_collision_design_cs_sweep as DS
    import collision_designs as CD
    from HAT_collision_on_off_position import REBUILT_HATCH as RB
    from site_layer import hat_topo_version as tv
    lv = pd.read_csv(tv.product_dir(START) / "survey-crests" / VOLUMES)
    st = pd.read_csv(tv.product_dir(START) / "survey-crests" / "survey_dune_crests.csv",
                     usecols=["survey", "model_state"]).drop_duplicates("survey")
    state = {int(r.survey): pd.Timestamp(r.model_state).year for r in st.itertuples()}
    base = lv[lv.survey == min(state)].set_index(["gis", "cell"]).volume_m3_per_m
    lv = lv[lv.survey > min(state)].copy()
    lv["obs"] = lv.volume_m3_per_m.to_numpy() - base.reindex(pd.MultiIndex.from_frame(lv[["gis", "cell"]])).to_numpy()
    lv = lv[lv.gis.between(GIS_LO, GIS_HI) & lv.obs.notna()
            & ~(lv.gis.between(*RB["gis"]) & lv.survey.isin(RB["surveys"]))]
    years = {y: s for y, s in state.items() if y in set(lv.survey)}

    def vol_score(v):
        dv = model_volume_change(v, years)
        lv["mod"] = [dv[s][g - 1, c] for s, g, c in zip(lv.survey, lv.gis, lv.cell)]
        dom = lv.groupby(["survey", "gis"])[["obs", "mod"]].median()
        e = dom["mod"] - dom.obs
        return float(np.sqrt((e ** 2).mean())), float(e.mean())

    rows = []
    for ref in (d for d in CD.DESIGNS if d.erosion == "nocollision"):
        m = DS.score(ref.run, 0.0)
        r, b = vol_score(ref.run)
        rows.append(dict(design=ref.name, erosion=ref.erosion, start=ref.start, cs=0.0, run=ref.run,
                         crest_rmse_m=m["crest_rmse_m"], volume_rmse_m3_m=r, volume_bias_m3_m=b))
        print(f"{ref.name}: crest {m['crest_rmse_m']:.3f} m, volume {r:.1f} m3/m (bias {b:+.1f})", flush=True)
    for d in (d for d in CD.DESIGNS if d.erosion != "nocollision"):
        k = f"{d.erosion}_{d.start}"
        for c in DS.H.sweep_cs(k):
            v = DS.H.DESIGN_SWEEP[(k, c)]
            if not DS.R._done(v):
                continue
            m = DS.score(v, c * 1e-6)
            r, b = vol_score(v)
            rows.append(dict(design=", ".join(d.name.split(", ")[:2]), erosion=d.erosion, start=d.start,
                             cs=c * 1e-6, run=v, crest_rmse_m=m["crest_rmse_m"], volume_rmse_m3_m=r,
                             volume_bias_m3_m=b))
            print(f"{d.name.split(', Cs')[0]}, Cs {c}e-6: crest {m['crest_rmse_m']:.3f} m, volume {r:.1f} m3/m "
                  f"(bias {b:+.1f})", flush=True)
    t = pd.DataFrame(rows)
    ref = t[t.erosion == "nocollision"].set_index("start")
    t["test1_beats_no_collision"] = t.crest_rmse_m < t.start.map(ref.crest_rmse_m)
    t["next_cs_crest_change_m"] = t.groupby(["erosion", "start"]).crest_rmse_m.transform(lambda x: x.shift(-1) - x)
    t["test2_stable"] = t.next_cs_crest_change_m.abs() < STABLE_M
    t["test3_volume_beats_no_collision"] = t.volume_rmse_m3_m < t.start.map(ref.volume_rmse_m3_m)
    t["passes_all"] = t.test1_beats_no_collision & t.test2_stable & t.test3_volume_beats_no_collision
    tests = ["test1_beats_no_collision", "test2_stable", "test3_volume_beats_no_collision", "passes_all"]
    t[tests] = t[tests].astype(object)
    t.loc[t.erosion == "nocollision", tests] = None
    t["cs"] = t.cs.map(lambda x: f"{x:.0e}" if x else "-")
    out = CD.ALL_DIR / "supporting" / OUT
    t.round(4).to_csv(out, index=False)
    print(t[["design", "cs", "crest_rmse_m", "volume_rmse_m3_m", "volume_bias_m3_m", "test1_beats_no_collision",
             "test2_stable", "test3_volume_beats_no_collision", "passes_all"]].round(3).to_string(index=False))
    print(out)


HALF_RUNS = (("nocollision_equal", None), ("lowerboth_equal", 2), ("faceretreat_equal", 3),
             ("faceretreat_equal", 10), ("lowerfront_equal", 1), ("lowerfront_equal", 2),
             ("lowerfront_equal", 3), ("lowerfront_equal", 30), ("lowerfront_equal", 200))
HALF_ZONES = {"Island (GIS 6-90)": range(6, 91), "Buxton (9-14)": range(9, 15), "Avon south (23-26)": range(23, 27),
              "Rodanthe (82-89)": range(82, 90), "Rest (no hotspot)": [g for g in range(6, 91) if not (
                  9 <= g <= 14 or 23 <= g <= 26 or 82 <= g <= 89)]}
HALVES_OUT = "seaward_vs_landward_1996_2009.csv"


def halves():
    """Did the real dune lose its seaward half and keep its landward half (as Lower front row assumes), 1996 to
    2009? Lidar: the change in mean height above the berm of each half of the dune window, 1996 survey to 2009
    survey. Model: the change in the front and back row, start to the state the 2009 survey is compared with.
    Per GIS domain the median over cells both have, then the mean over the domains of each zone."""
    import collision_designs as CD
    import HAT_collision_design_cs_sweep as DS
    from site_layer import hat_topo_version as tv
    from cascade_pipeline.saved_runs import load_saved_cascade
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    lv = pd.read_csv(tv.product_dir(START) / "survey-crests" / VOLUMES).set_index(["survey", "gis", "cell"])
    st = pd.read_csv(tv.product_dir(START) / "survey-crests" / "survey_dune_crests.csv",
                     usecols=["survey", "model_state"]).drop_duplicates("survey")
    end_state = pd.Timestamp(st[st.survey == 2009].model_state.iloc[0]).year
    obs = (lv.loc[2009][["seaward_half_m", "landward_half_m"]] - lv.loc[1996][["seaward_half_m", "landward_half_m"]])
    obs = obs.dropna()
    keep = obs.index
    rows = [dict(source="Lidar", zone=z, seaward_or_front_m=float(obs.loc[list(gs)].groupby("gis").seaward_half_m
                                                                    .median().mean()),
                 landward_or_back_m=float(obs.loc[list(gs)].groupby("gis").landward_half_m.median().mean()))
            for z, gs in HALF_ZONES.items()]
    for k, c in HALF_RUNS:
        v = next(d.run for d in CD.DESIGNS if f"{d.erosion}_{d.start}" == k) if c is None else DS.H.DESIGN_SWEEP[(k, c)]
        rd = DS.H.run_dir(v)
        first = int(json.loads(next(rd.glob("*_run_metadata.json")).read_text(encoding="utf-8"))
                    ["period"]["first_model_year"])
        cas = load_saved_cascade(rd / f"{rd.name}.npz")
        dm = DOM.for_total(len(cas.barrier3d))
        recs = []
        for g in range(6, 91):
            dd = np.asarray(cas.barrier3d[dm.gis_to_pad(g)].DuneDomain, float) * 10
            ch = dd[end_state - first] - dd[0]                                        # (50, 2) m
            cells = [cc for gg, cc in keep if gg == g]
            if cells:
                recs.append(dict(gis=g, front=float(np.median(ch[cells, 0])), back=float(np.median(ch[cells, 1]))))
        m = pd.DataFrame(recs).set_index("gis")
        lab = CD.EROSION[k.split("_")[0]][0] + ("" if c is None else f", Cs {DS.cs_label(c)}")
        for z, gs in HALF_ZONES.items():
            mm = m.loc[[g for g in gs if g in m.index]]
            rows.append(dict(source=lab, zone=z, seaward_or_front_m=float(mm.front.mean()),
                             landward_or_back_m=float(mm.back.mean())))
        print(f"loaded {lab}", flush=True)
    t = pd.DataFrame(rows)
    t["seaward_minus_landward_m"] = t.seaward_or_front_m - t.landward_or_back_m
    out = CD.ALL_DIR / "supporting" / HALVES_OUT
    t.round(3).to_csv(out, index=False)
    for z in HALF_ZONES:
        print(f"\n{z}: change 1996-2009, m (seaward half / front row, landward half / back row)")
        print(t[t.zone == z][["source", "seaward_or_front_m", "landward_or_back_m", "seaward_minus_landward_m"]]
              .round(2).to_string(index=False))
    print(out)


if __name__ == "__main__":
    {"lidar": lidar, "score": score, "halves": halves}[sys.argv[1] if len(sys.argv) > 1 else "score"]()
