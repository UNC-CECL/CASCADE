"""
Does the lidar show the sand the model overwashes? (2026-11-19) Elevation change in the strip behind the dune between
two surveys, per GIS domain, against the run's overwash over the same model years: an independent check on the
"false alarms" of the image comparison, which may be thin, cleared or reworked washover the photos miss.

    python scripts/real-world-cases/Hatteras/dune-collision/4-hindcast/3-on-vs-off/HAT_overwash_lidar_deposition.py

Lidar: per 10 m alongshore cell, the mean change in ground height over the BEHIND_M strip landward of the survey's
dune window (the landward edge of the window the crests are picked in, survey_window_shifts.csv; the later survey's
window), times BEHIND_M, as m3 per m of shore; per domain the median over its cells. Water and unsurveyed cells are
left out. Model: the run's overwash (QowTS, m3/m) summed over the model years between the two surveys' model
states. Pairs: 2000 -> 2004 (holds Isabel, the control) and 2005 -> 2009 (the period of the unconfirmed overwash).
Each survey's datum check: the median change over all land cells of the island strip, printed.

Writes 3-collision-design/0-all-designs/overwash_lidar_deposition_<pair>.png and supporting/...csv.

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
for p in (_HAT, _HERE.parent, _HAT / "input_prep" / "1-barrier3d-domains" / "2-survey-crests"):
    sys.path.insert(0, str(p))
import collision_designs as CD  # noqa: E402
import HAT_collision_row_mode_sweep as R  # noqa: E402

H = R.H

# --- CONFIG ------------------------------------------------------------------
START = "1996-start"
PAIRS = ((2000, 2004), (2005, 2009))
BEHIND_M = 150                       # strip landward of the dune window, m
GIS_LO, GIS_HI = 1, 90
RUNS = (("No collision, equal rows", "cal_step2_growth_noscrape_vdefault", CD.EROSION["nocollision"][2]),
        ("Face retreat, equal rows, Cs 3e-6", H.DESIGN_SWEEP[("faceretreat_equal", 3)], CD.EROSION["faceretreat"][2]),
        ("Lower front row, equal rows, Cs 2e-4", H.DESIGN_SWEEP[("lowerfront_equal", 200)], CD.EROSION["lowerfront"][2]))
MODEL_OVERWASH_M3 = 50.0             # model overwash counted as "substantial" (a visible fan), m3/m over the pair
LIDAR_GAIN_M3 = 15.0                 # lidar gain behind the dune counted as deposition, m3/m (0.1 m over 150 m)
# -----------------------------------------------------------------------------


def lidar_strips(years):
    """{year: {gis: (z behind the dune, (50, cells) m NAVD88, NaN where water/unsurveyed)}} in the later survey's
    windows, plus each survey's land cells for the datum check."""
    import rasterio
    import HAT_survey_dune_crests as SC
    X, tv = SC.X, SC.tv
    prod, surveys = SC.STARTS[START]
    windows = json.loads(tv.picks_dir(START).joinpath(
        f"HAT_dune_search_windows_{tv.topo_dirs(START)[2]}.json").read_text(encoding="utf-8"))
    frames = {g: SC.frame(START, g, windows) for g in range(GIS_LO, GIS_HI + 1)}
    grids = {g: SC.clip_grid(SC.product(prod).gapfill_1m / f"clip_domain_{g}_filled.tif") for g in range(1, 91)}
    out = {}
    for year in years:
        kind, path = surveys[year][:2]
        files = sorted(path.parent.glob(path.name)) if "*" in path.name else [path]
        srcs = [rasterio.open(f) for f in files]
        try:
            unit = "m" if kind == "mosaic" else SC.unit_of(srcs, START, grids)[0]
            scale = SC.US_FT if unit == "ft" else 1.0
            shift = SC.SURVEY_VERTICAL_OFFSET_M.get(START, {}).get(year, 0.0)
            out[year] = {}
            for g in range(GIS_LO, GIS_HI + 1):
                shear, c0, *_ = frames[g]
                on = SC.mosaic_on_clip(srcs[0], grids[g]) if kind == "mosaic" else SC.raster_on_clip(srcs, grids[g])
                z = SC.profiles(SC.downsample(on * scale + shift), shear, c0)[:50]
                z = np.where(z <= X.SENTINEL_WATER_M + 1e-9, np.nan, z) + X.MHW_M        # water, nodata -> NaN
                out[year][g] = z
        finally:
            for s in srcs:
                s.close()
        print(f"  {year} read", flush=True)
    return out, frames, SC


def behind_change(z0, z1, i1s, n):
    """Per alongshore cell: mean change over the n cells landward of the window edge i1s, m."""
    d = np.full(50, np.nan)
    for c in range(min(50, z0.shape[0], z1.shape[0])):
        a = int(i1s[c])
        b = min(a + n, z0.shape[1], z1.shape[1])
        if b > a:
            dz = z1[c, a:b] - z0[c, a:b]
            if np.isfinite(dz).sum() >= n // 2:
                d[c] = np.nanmean(dz)
    return d


def model_overwash(v, y0, y1):
    """Per GIS domain: the run's overwash summed over the model years after survey y0's model state up to survey y1's
    (QowTS rows: row 0 the initial state, row k model year first + k - 1), m3/m."""
    rd = H.run_dir(v)
    meta = json.loads(next(rd.glob("*_run_metadata.json")).read_text(encoding="utf-8"))
    first = int(meta["period"]["first_model_year"])
    q = np.load(rd / f"{rd.name}_overwash.npy", allow_pickle=True).astype(float)
    nbuf = (q.shape[1] - 90) // 2
    real = q[:, nbuf: nbuf + 90]
    st = pd.read_csv(H.SURVEY_CRESTS, usecols=["survey", "model_state"]).drop_duplicates("survey").set_index("survey")
    s0, s1 = pd.Timestamp(st.model_state[y0]).year, pd.Timestamp(st.model_state[y1]).year   # 1 Jan states
    rows = [k for k in range(1, real.shape[0]) if s0 <= first + k - 1 < s1]
    return real[rows].sum(axis=0), (first + rows[0] - 1, first + rows[-1] - 1)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer import hat_figure_style as S
    S.apply_style()
    years = sorted({y for p in PAIRS for y in p})
    z, frames, SC = lidar_strips(years)
    n = BEHIND_M // 10
    for y0, y1 in PAIRS:
        rows = []
        for g in range(GIS_LO, GIS_HI + 1):
            _, _, i0, i1 = frames[g]
            _, i1s = SC.survey_windows(START, y1, g, i0, i1)
            d = behind_change(z[y0][g], z[y1][g], i1s, n)
            rows.append(dict(gis=g, lidar_change_m=float(np.nanmedian(d)) if np.isfinite(d).any() else np.nan))
        t = pd.DataFrame(rows).set_index("gis")
        t["lidar_m3_per_m"] = t.lidar_change_m * BEHIND_M
        land = np.concatenate([(z[y1][g] - z[y0][g]).ravel() for g in range(GIS_LO, GIS_HI + 1)])
        print(f"\n{y0} -> {y1}: datum check, median change over all land cells {np.nanmedian(land):+.2f} m")
        span = None
        for lab, v, _ in RUNS:
            ow, span = model_overwash(v, y0, y1)
            t[lab] = ow
        out = CD.ALL_DIR / "supporting" / f"overwash_lidar_deposition_{y0}_{y1}.csv"
        t.round(2).to_csv(out)
        gain = t.lidar_m3_per_m > LIDAR_GAIN_M3
        print(f"model years {span[0]}-{span[1]}; lidar gain > {LIDAR_GAIN_M3:g} m3/m behind the dune in "
              f"{int(gain.sum())} of {int(t.lidar_m3_per_m.notna().sum())} domains")
        for lab, v, _ in RUNS:
            sub = t[lab] > MODEL_OVERWASH_M3
            both = (sub & gain).sum()
            r = t[[lab, "lidar_m3_per_m"]].dropna().corr().iloc[0, 1]
            print(f"  {lab:38s}: substantial overwash in {int(sub.sum()):2d} domains; lidar gain in {int(both):2d} "
                  f"of them ({100 * both / max(sub.sum(), 1):.0f}%); lidar gain where it had none: "
                  f"{int((gain & ~sub).sum())}; correlation {r:+.2f}")
        # figure: lidar change behind the dune, and each run's overwash, along the island
        fig, axs = plt.subplots(2, 1, figsize=(S.FIG_W_DOUBLE, 5.0), sharex=True, layout="constrained",
                                gridspec_kw={"height_ratios": [1, 1]})
        ax = axs[0]
        col = np.where(t.lidar_m3_per_m >= 0, "#2C7FB8", "#D7301F")
        ax.bar(t.index, t.lidar_m3_per_m, color=col, width=0.85)
        ax.axhline(0, color="0.3", lw=0.7)
        ax.axhline(LIDAR_GAIN_M3, color="0.5", lw=0.6, ls=":")
        ax.set_ylabel("Lidar change behind\nthe dune (m³/m)", fontsize=8.5)
        ax.set_title(f"Lidar, {y0} → {y1} survey: ground change in the {BEHIND_M} m behind the dune "
                     "(blue gain, red loss)", fontsize=9, loc="left")
        ax2 = axs[1]
        w = 0.8 / len(RUNS)
        for k, (lab, v, c) in enumerate(RUNS):
            ax2.bar(t.index + (k - 1) * w, t[lab], width=w, color=c, label=lab)
        ax2.set_ylabel("Model overwash\n(m³/m)", fontsize=8.5)
        ax2.set_title(f"Model overwash, model years {span[0]}–{span[1]}", fontsize=9, loc="left")
        ax2.legend(fontsize=7.5, frameon=False, ncol=3, loc="upper center")
        ax2.set_xlabel("GIS domain (south → north)")
        for a in axs:
            a.grid(axis="y", color="0.92", lw=0.5)
            a.spines[["top", "right"]].set_visible(False)
            a.set_xlim(0.3, 90.7)
        o = CD.ALL_DIR / f"overwash_lidar_deposition_{y0}_{y1}.png"
        S.save(fig, o, bbox_inches="tight", close=True)
        S.record_caption(o, (
            f"Top: the lidar change in ground height between the {y0} and {y1} surveys over the {BEHIND_M} m "
            "landward of each profile's dune window (the window the crests are picked in, the later survey's), as m3 "
            f"per m of shore, the median over each GIS domain's 10 m cells; dotted: {LIDAR_GAIN_M3:g} m3/m. Bottom: "
            f"the model's overwash over model years {span[0]}-{span[1]} for " + "; ".join(
                f"{lab} ({v})" for lab, v, _ in RUNS) + ". Lidar change behind the dune also includes vegetation "
            "growth and loss, dune growth at the window edge, road clearing and survey differences."))
        print(f"wrote {o}")


if __name__ == "__main__":
    main()
