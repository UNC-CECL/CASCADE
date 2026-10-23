"""
The seven calibration domains followed into the 2009-2025 test period, one folder per domain.

    python scripts/hatteras_ms/experiments/HAT_test_period_domain_figures.py

The 2009 counterparts of the 1996 start's per-domain observed-dune figures
(HAT_stable_dune_figures.py, 2-observed-dunes/gis<N>/): the same seven domains
(GIS 20-22 and 59-62), from the 2009 lidar to the 2019 lidar, with 2014-2018
between. Observed only; no model runs. Written to
data/.../2009-start/dune-growth-rate/figures/observed-dunes/gis<N>/:
  test_dune_gis<N>_dem          the 2009 and 2019 1 m DEMs and their difference, north up, with the 2009 dune
                                window and each row's picked crest
  test_dune_gis<N>_profile      2009 and 2019 cross-shore profiles
  test_dune_gis<N>_surveys      the crest through all six surveys, 2009-2019
  test_dune_gis<N>_alongshore   the crest along the domain in every survey
  test_dune_gis<N>_crest_shift  how far the picked crest moved cross-shore, 2009 to 2019
Surveys are cut exactly as the survey-crest table (HAT_survey_dune_crests.py): the 1 m mosaic on the 2009-2014
product's clip grid, shifted to the 2009 lidar on NC-12, in the 2009 start's frame and dune window. The script
stops unless its 2014-2019 crests reproduce that table.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-21
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import rasterio  # noqa: E402

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_stable_dune_figures as F  # noqa: E402  (the 1996 per-domain figures: frame, colours, helpers)
import HAT_test_period_growth_rate as T  # noqa: E402
from site_layer import hat_figure_style as fs  # noqa: E402
from site_layer.hat_elevation_products import product  # noqa: E402

S, X, tv = F.S, F.X, F.tv

# --- CONFIG ------------------------------------------------------------------
OUT = T.OUT_DIR / "figures" / "observed-dunes"
START = "2009-start"
PROD, SURVEYS = S.STARTS[START]             # "2009-2014", {2014: (kind, path, ...), ...}
SHIFT = S.SURVEY_VERTICAL_OFFSET_M[START]   # each survey moved to the 2009 lidar on NC-12
FIRST, LAST = 2009, 2019
C_FIRST, C_LAST = F.C09, fs.C_1997          # house colours: 2009 amber, the test period blue
STORMS = {"Irene": 2011.65, "Sandy": 2012.82, "Matthew": 2016.77, "Florence": 2018.71, "Dorian": 2019.68}
# -----------------------------------------------------------------------------


def save(fig, name, caption):
    path = OUT / f"gis{name.split('_')[2][3:]}" / f"{name}.png"
    fs.save(fig, path, close=True)
    fs.record_caption(path, caption + " Producer: HAT_test_period_domain_figures.py.")
    print(path, flush=True)


# --- data ----------------------------------------------------------------------------
_WINDOWS = {}


def frame(gis):
    if not _WINDOWS:
        _WINDOWS.update(json.loads(tv.picks_dir(START).joinpath(
            f"HAT_dune_search_windows_{tv.topo_dirs(START)[2]}.json").read_text(encoding="utf-8")))
    return S.frame(START, gis, _WINDOWS)


# One survey's 1 m DEM (m NAVD88) on the clip grid; NaN where that survey has no data. 2009 is the start's own
# lidar (not its 2014 gap fill), the later ones the mosaics shifted as in the survey-crest table
def dem_1m(year, gis):
    root = product(PROD).gapfill_1m
    with rasterio.open(root / f"clip_domain_{gis}_filled.tif") as s:
        z, t, grid = S.clean(s.read(1), s.nodata), s.transform, (s.transform, s.shape, s.crs)
    if year == FIRST:
        with rasterio.open(root / f"clip_domain_{gis}_survey.tif") as s:
            z[s.read(1) != FIRST] = np.nan
        return z, t
    with rasterio.open(SURVEYS[year][1]) as src:
        return S.mosaic_on_clip(src, grid) + SHIFT.get(year, 0.0), t


_PICKS = {}


# The crest picked in each 10 m row and its frame column (-1 where none), m NAVD88; as F.crest_picks. Also
# where the crest is cut by the window: on its last cell with higher ground just landward (the dune has moved
# out of the 2009 window, so the crest reads low). Most last-cell crests are real peaks and are kept
def crest_picks(year, gis):
    if (year, gis) not in _PICKS:
        shear, c0, i0, i1 = frame(gis)
        z = S.profiles(S.downsample(dem_1m(year, gis)[0]), shear, c0)
        height, edge, status = S.crests(z, i0, i1)
        above = z > X.BEACH_START_THR_M
        sb = np.where(above.any(axis=1), above.argmax(axis=1), -1)
        pos = np.full(z.shape[0], -1)
        for i in np.flatnonzero(status == "ok"):
            a = max(i0, int(sb[i])) if (X.CLIP_WINDOW_TO_BEACH and sb[i] >= 0) else i0
            win = z[i, a:min(i1, z.shape[1])]
            pos[i] = a + int(np.argmax(np.where(win > X.SENTINEL_WATER_M + 1e-9, win, -np.inf)))
        nxt = z[:, i1] if i1 < z.shape[1] else np.full(z.shape[0], -np.inf)
        _PICKS[(year, gis)] = (pos, height + X.MHW_M, edge & (nxt > height))
    return _PICKS[(year, gis)]


# Crest per cell and survey, 2009-2019, checked against the survey-crest table; and the rows cut by the window
def cell_crests(gis):
    t = pd.read_csv(tv.product_dir(START) / "survey-crests" / "survey_dune_crests.csv")
    t = t[t.gis == gis]
    w = t.pivot_table(index="cell", columns="survey", values="crest_m_navd88")
    for year in w.columns:                    # this script's cut must be the table's
        h = crest_picks(year, gis)[1]
        d = np.abs(w[year].to_numpy() - h[w.index])
        if np.nanmax(d) > 0.002 or (np.isnan(w[year].to_numpy()) != np.isnan(h[w.index])).any():
            raise SystemExit(f"GIS {gis} {year}: crests differ from survey_dune_crests.csv by {np.nanmax(d):.3f} m")
    w.insert(0, FIRST, crest_picks(FIRST, gis)[1][w.index])
    cut = pd.DataFrame({y: crest_picks(y, gis)[2][w.index] for y in w.columns}, index=w.index)
    return w, cut


def cut_note(cut):
    n = cut.sum()
    n = n[n > 0]
    if n.empty:
        return ""
    return (" Crests cut by the window (on its landward edge with higher ground just landward, so they read "
            "low): " + ", ".join(f"{c} of 50 rows in {y}" for y, c in n.items()) + ".")


def profile_array(year, gis):
    shear, c0, i0, i1 = frame(gis)
    return F.to_frame(S.downsample(dem_1m(year, gis)[0]), shear, c0), (i0, i1)


def stretch_of(gis):
    return F.stretch_of(gis)


# --- figures -------------------------------------------------------------------------
def fig_dem(gis):
    za, t = dem_1m(FIRST, gis)
    zb, _ = dem_1m(LAST, gis)
    pa, _, _ = crest_picks(FIRST, gis)
    pb, _, _ = crest_picks(LAST, gis)
    shear, c0, i0, i1 = frame(gis)
    rows = np.arange(len(pa))
    xs_lo = [F.frame_to_xy(i, i1 - 1, shear, c0, t)[0] for i in rows]
    xs_hi = [F.frame_to_xy(i, i0, shear, c0, t)[0] for i in rows]
    x0, x1 = min(xs_lo) - 80, max(xs_hi) + 40
    c_lo, c_hi = int(x0 - t.c), int(x1 - t.c)
    ext = (t.c + c_lo, t.c + c_hi, t.f - za.shape[0], t.f)
    cmap, _, _ = fs.elevation_cmap()
    norm = matplotlib.colors.Normalize(0.0, 8.0)
    ls = matplotlib.colors.LightSource(azdeg=315, altdeg=45)
    fig, axes = plt.subplots(1, 3, figsize=fs.figsize("double", height=5.2), layout="constrained",
                             sharex=True, sharey=True)
    titles = (f"Observed: {FIRST} lidar (USACE)", f"Observed: {LAST} lidar (DUNEX)", f"{LAST} minus {FIRST}")
    for k, (ax, z, title) in enumerate(zip(axes, (za, zb, zb - za), titles)):
        ax.set_facecolor("0.85")
        if k < 2:
            zc = z[:, c_lo:c_hi] - X.MHW_M
            ax.imshow(ls.hillshade(np.nan_to_num(zc, nan=np.nanmin(zc)), vert_exag=2), cmap="gray",
                      extent=ext, interpolation="none")
            im = ax.imshow(np.ma.masked_invalid(zc), cmap=cmap, norm=norm, alpha=0.75, extent=ext,
                           interpolation="none")
        else:
            imd = ax.imshow(np.ma.masked_invalid(z[:, c_lo:c_hi]), cmap="RdBu", vmin=-2, vmax=2, extent=ext,
                            interpolation="none")
        for j, side in ((i0, 1), (i1 - 1, -1)):
            xy = np.array([F.frame_to_xy(i, j, shear, c0, t) for i in rows])
            ax.step(xy[:, 0] + side * 5, xy[:, 1], where="mid", color=fs.INK, lw=0.7, ls="--")
        for pos, col, year in ((pa, C_FIRST, FIRST), (pb, C_LAST, LAST)):
            if (k == 0 and year == LAST) or (k == 1 and year == FIRST):
                continue
            ok = pos >= 0
            xy = np.array([F.frame_to_xy(i, pos[i], shear, c0, t) for i in rows[ok]]).reshape(-1, 2)
            ax.plot(xy[:, 0], xy[:, 1], ls="none", marker="o", ms=2.6, mfc=col, mec="white", mew=0.4,
                    label=f"Crest picked from the {year} lidar")
        ax.set_title(title, fontsize=7)
        ax.set_aspect("equal")
        ax.set(xlim=(ext[0], ext[1]), ylim=(ext[2], ext[3]), xticks=[], yticks=[])
        fs._letter_inside(ax, k)
    a = axes[0]
    a.plot([ext[0] + 10, ext[0] + 60], [ext[2] + 15] * 2, color=fs.INK, lw=2)
    a.text(ext[0] + 35, ext[2] + 22, "50 m", ha="center", fontsize=6)
    a.annotate("", (ext[0] + 20, ext[3] - 40), xytext=(ext[0] + 20, ext[3] - 80),
               arrowprops=dict(arrowstyle="-|>", color=fs.INK, lw=0.8))
    a.text(ext[0] + 20, ext[3] - 36, "N", ha="center", va="bottom", fontsize=7)
    a.text(ext[1] - 4, ext[2] + 15, "ocean →", ha="right", fontsize=6.5, color=fs.INK,
           bbox=dict(facecolor="white", edgecolor="none", pad=1, alpha=0.8))
    fig.colorbar(im, ax=axes[:2], orientation="horizontal", shrink=0.6, label="Elevation (m MHW)")
    fig.colorbar(imd, ax=axes[2], orientation="horizontal", shrink=0.9, label="Change (m)")
    handles = [h for ax in axes[:2] for h in ax.get_legend_handles_labels()[0]]
    labels = [lab for ax in axes[:2] for lab in ax.get_legend_handles_labels()[1]]
    handles.append(matplotlib.lines.Line2D([], [], color=fs.INK, lw=0.7, ls="--"))
    labels.append("Model dune window (2009 start)")
    fig.legend(handles, labels, frameon=False, fontsize=6.5, loc="outside lower center", ncol=3)
    fig.suptitle(f"GIS {gis} ({stretch_of(gis)})", fontsize=8, fontweight="bold", x=0.01, ha="left")
    both = (pa >= 0) & (pb >= 0)
    shift = (pb[both] - pa[both]) * fs.CELL_M
    gone = int(((pa >= 0) & (pb < 0)).sum())
    save(fig, f"test_dune_gis{gis}_dem",
         f"The same dune in the 2009 and 2019 surveys, GIS {gis} (in the {stretch_of(gis)} stretch). (a) The 2009 "
         f"lidar (USACE), the model's starting surface, and (b) the 2019 lidar (DUNEX, lowered 0.12 m to the 2009 "
         f"lidar on NC-12) at 1 m, north up, true scale, over a hillshade; grey where a survey has no data. "
         f"(c) 2019 minus 2009: blue gain, red loss. Dashed: the 2009 start's dune window, the band of 10 m cells "
         f"each row's crest is searched in. Dots: the crest picked in each 10 m alongshore row from that year's "
         f"DEM, averaged to 10 m as the model sees it. {np.mean(np.abs(shift) <= fs.CELL_M):.0%} of the "
         f"{both.sum()} rows picked in both years are within one 10 m cell of each other"
         + (f"; {gone} rows have no 2019 crest (the dune window is now beach)." if gone else ".")
         + cut_note(pd.DataFrame({y: crest_picks(y, gis)[2] for y in (FIRST, LAST)})))


def fig_profile(gis):
    fig, ax = plt.subplots(figsize=fs.figsize("single", height=2.8), layout="constrained")
    xs = None
    for year, col in ((FIRST, C_FIRST), (LAST, C_LAST)):
        z, (i0, i1) = profile_array(year, gis)
        med = np.nanmedian(z, axis=0)
        lo, hi = np.nanpercentile(z, 25, axis=0), np.nanpercentile(z, 75, axis=0)
        if xs is None:                                       # x = 0 at the 2009 shoreline (MHW crossing)
            j0 = int(np.argmax(med > X.MHW_M))
            xs = (np.arange(z.shape[1]) - j0) * fs.CELL_M
            win = ((i0 - j0) * fs.CELL_M, (i1 - j0) * fs.CELL_M)
        ax.fill_between(xs, lo, hi, color=col, alpha=0.18, lw=0)
        ax.plot(xs, med, color=col, lw=1.6, label=f"Observed: {year} lidar DEM")
    ax.axvspan(*win, color="0.93", zorder=0)
    ax.axhline(X.MHW_M, color=fs.INK_MUTED, lw=0.6, ls=":")
    ax.text(win[0] + 3, -0.6, "dune window", fontsize=6, color=fs.INK_MUTED)
    ax.set(xlim=(-40, win[1] + 60), ylim=(-1, 9.5), xlabel=f"Distance landward of the {FIRST} shoreline (m)",
           ylabel="Elevation (m NAVD88)")
    ax.text(0.97, 0.95, f"GIS {gis}", transform=ax.transAxes, ha="right", va="top", fontsize=7, fontweight="bold")
    ax.legend(frameon=False, fontsize=6.5, loc="upper left")
    F.finish(ax)
    save(fig, f"test_dune_gis{gis}_profile",
         f"Cross-shore profile at GIS {gis} (in the {stretch_of(gis)} stretch) from the 2009 and 2019 lidar DEMs, "
         f"in the 2009 start's frame (10 m cells; line: median of the domain's 50 alongshore rows, band: 25th-75th "
         f"percentile; only cells each survey covers). Grey: the dune window the model reads its dune from.")


def fig_surveys(gis):
    w, cut = cell_crests(gis)
    w = w.dropna()
    years = [int(y) for y in w.columns]
    lab = stretch_of(gis)
    fig, ax = plt.subplots(figsize=fs.figsize("single", height=2.8), layout="constrained")
    ax.fill_between(years, w.quantile(0.25).to_numpy(), w.quantile(0.75).to_numpy(), color=F.STABLE_C[lab],
                    alpha=0.2, lw=0, label="Observed: middle half of the cells")
    ax.plot(years, w.median().to_numpy(), color=F.STABLE_C[lab], lw=2, marker="o", ms=4,
            label=f"Observed: GIS {gis}, median")
    for name, x in STORMS.items():
        ax.axvline(x, color=fs.INK_MUTED, lw=0.6, ls="--")
        ax.text(x + 0.08, 9.25, name, fontsize=6, color=fs.INK_MUTED, rotation=90, va="top")
    ax.set(xlabel="Survey year", ylabel="Dune crest (m NAVD88)", xlim=(2008.5, 2020.0), ylim=(4, 9.5))
    ax.set_xticks(years)
    ax.legend(frameon=False, fontsize=6.5, loc="lower right")
    F.finish(ax)
    rise = w[LAST].median() - w[FIRST].median()
    save(fig, f"test_dune_gis{gis}_surveys",
         f"Dune crest at GIS {gis} (in the {lab} stretch) through the six lidar surveys, 2009-2019: median and "
         f"middle half of the {len(w)} 10 m cells with a dune in every survey. The median changed {rise:+.2f} m "
         f"from 2009 to 2019. All are 1 m DEMs, 2014-2019 shifted to the "
         f"2009 lidar on NC-12. Dashed: storms between the surveys; the 2019 survey was flown in September-October, "
         f"just after Dorian." + cut_note(cut))


def fig_alongshore(gis):
    lab = stretch_of(gis)
    w, cut = cell_crests(gis)
    x = w.index * fs.CELL_M + 5
    fig, ax = plt.subplots(figsize=fs.figsize("double", height=2.5), layout="constrained")
    for year in w.columns:
        if year not in (FIRST, LAST):
            ax.plot(x, w[year], color="0.8", lw=0.6, zorder=1)
    for year, col in ((FIRST, C_FIRST), (LAST, C_LAST)):
        ax.plot(x, w[year], color=col, lw=1.4, label=f"Observed: {year} lidar DEM (1 m)", zorder=3)
    ax.plot([], [], color="0.8", lw=0.6, label="Observed: the 2014-2018 surveys (1 m)")
    both = w[[FIRST, LAST]].dropna()
    r = np.corrcoef(both[FIRST], both[LAST])[0, 1]
    rise = np.median(both[LAST] - both[FIRST])
    ax.text(0.99, 0.97, f"GIS {gis}: {FIRST} vs {LAST} alongshore r = {r:.2f}, median change {rise:+.2f} m",
            transform=ax.transAxes, ha="right", va="top", fontsize=6.5, color=F.STABLE_C[lab],
            fontweight="bold", bbox=dict(facecolor="white", edgecolor="none", pad=1, alpha=0.8))
    ax.set(xlim=(0, 500), ylim=(3, 9.8), ylabel="Dune crest (m NAVD88)",
           xlabel="Alongshore distance, south to north (m)")
    F.finish(ax)
    fig.legend(*ax.get_legend_handles_labels(), frameon=False, fontsize=6.5, loc="outside lower center", ncol=3)
    gone = int(w[LAST].isna().sum() - w[FIRST].isna().sum())
    save(fig, f"test_dune_gis{gis}_alongshore",
         f"Dune crest along GIS {gis} (in the {lab} stretch), one value per 10 m model row (the highest 10 m cell "
         f"in the 2009 start's dune window). Amber and blue: the 2009 and 2019 1 m DEMs; grey: 2014, 2016, 2017 "
         f"and 2018. 2009 and 2019 correlate alongshore at r = {r:.2f} and the median row changed {rise:+.2f} m."
         + (f" Gaps: {gone} rows where the 2009 dune position is beach by 2019." if gone > 0 else "")
         + cut_note(cut))


def fig_crest_shift(gis):
    pa, _, _ = crest_picks(FIRST, gis)
    pb, _, _ = crest_picks(LAST, gis)
    ok = (pa >= 0) & (pb >= 0)
    d = (pb[ok] - pa[ok]) * fs.CELL_M
    lab = stretch_of(gis)
    fig, ax = plt.subplots(figsize=fs.figsize("single", height=2.6), layout="constrained")
    ax.hist(d, bins=np.arange(-65, 70, 10), color=F.STABLE_C[lab], rwidth=0.85, label=f"Observed: GIS {gis}")
    within = np.mean(np.abs(d) <= fs.CELL_M)
    ax.axvline(0, color=fs.INK, lw=0.6, ls="--")
    ax.text(0.98, 0.95, f"{within:.0%} of {len(d)} rows within\none 10 m cell", transform=ax.transAxes,
            ha="right", va="top", fontsize=6.5)
    ax.set(xlabel=f"Crest position, {LAST} minus {FIRST} (m, + landward)", ylabel="10 m rows")
    ax.legend(frameon=False, fontsize=6.5, loc="upper left")
    F.finish(ax)
    save(fig, f"test_dune_gis{gis}_crest_shift",
         f"How far the picked dune crest moved cross-shore between the 2009 and 2019 lidar at GIS {gis} (in the "
         f"{lab} stretch), per 10 m row picked in both years. 0 means the same 10 m cell; {within:.0%} of rows are "
         f"within one cell.")


def main():
    fs.apply_style()
    for gis in sorted(g for doms in F.STABLE.values() for g in doms):
        fig_dem(gis)
        fig_profile(gis)
        fig_surveys(gis)
        fig_alongshore(gis)
        fig_crest_shift(gis)


if __name__ == "__main__":
    main()
