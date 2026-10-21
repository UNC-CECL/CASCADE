"""
Why the stable dunes at GIS 59-62 and 20-22 can calibrate the dune growth rate r, one point per figure.

    python scripts/hatteras_ms/experiments/HAT_stable_dune_figures.py

Figures for the stable-dune study (HAT_stable_dune_growth_rate.py, run it first),
written to its figures/ folder (data/.../1996-start/dune-growth-rate/ since 2026-10-21), one subfolder per
topic: 1-selection (fig1, fig3), 2-observed-dunes (fig2, fig4, fig7-9), 3-model-fit (fig5, fig6). In reading order:
  stable_dune_fig1_where       shoreline change, dune change and overwash per domain, 1996-2009:
                               the two stretches meet all three conditions
  stable_dune_fig2_profiles    1996 and 2009 DEM profiles at the stable stretches: dune and beach in place
  stable_dune_fig3_excluded    the same at two excluded domains: a growing beach, a retreating shoreline
  stable_dune_fig4_surveys     crests through all nine 1996-2009 lidar surveys: no storm losses
  stable_dune_fig5_model       model crest for several r against the surveys
  stable_dune_fig6_rate        (drawn by HAT_stable_dune_growth_rate.py) model minus lidar against r
  stable_dune_fig7_dem_gis<N>  one per stable domain: the 1996 and 2009 1 m DEMs and their difference, north up,
                               with the dune window and each row's picked crest: the same dune in both
  stable_dune_fig8_alongshore  crest along the shore in every survey: 1996 and 2009 rise and dip together
  stable_dune_fig9_crest_shift how far the picked crest moved cross-shore, 1996 to 2009
    (fig7-9 alone, no model runs needed: HAT_stable_dune_figures.py dems)
Profiles: the 1 m gap-filled DEMs of the 1996 and 2009 starts (same clip grid),
averaged to 10 m as the starts are and put in the 1996 start's cross-shore frame
by the survey-crest code. A 10 m cell is drawn only where all its core 1 m
pixels come from that year's own survey (1996 ALACE, 2009 USACE); the 1996 start
is the 2009 DEM wherever ALACE has no data, so unmasked cells would fake stability.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-20
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
_REPO = next(p for p in _HERE.parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(_REPO / "scripts" / "input_prep" / "1-barrier3d-domains" / "2-survey-crests"))
import HAT_stable_dune_growth_rate as R  # noqa: E402  (points the launcher at the stable-dune study)
import HAT_survey_dune_crests as S  # noqa: E402
import HAT_dune_crest_vs_surveys as V  # noqa: E402
from site_layer import hat_figure_style as fs  # noqa: E402
from site_layer import hat_observed_rates as obsr  # noqa: E402
from site_layer import hat_overwash as ow  # noqa: E402
from site_layer.hat_elevation_products import product  # noqa: E402

X, tv = S.X, S.tv

# --- CONFIG ------------------------------------------------------------------
OUT = R.OUT_DIR / "figures"
# One subfolder per point the figures make, in reading order (2026-10-21)
TOPICS = {"stable_dune_fig1_where": "1-selection", "stable_dune_fig3_excluded": "1-selection",
          "stable_dune_fig5_model": "3-model-fit"}   # everything else: 2-observed-dunes
C96, C09 = fs.C_1984, "#c8880f"            # house period colours: 1996 red, 2009 amber
STABLE = R.GROUPS                           # {"GIS 59-62": (...), "GIS 20-22": (...)}
STABLE_C = R.GROUP_COLOURS
PROFILE_DOMAINS = (61, 21)                  # one per stable stretch
EXCLUDED = {44: "beach building out", 65: "shoreline retreating"}
SHORE_BAND = 0.5                            # m/yr: the stable-shoreline condition
MIN_CELLS = 25                              # domains with fewer fully surveyed cells are left blank in fig1
MODEL_RATES = {(0.0, 0.0): ("r = 0 (no growth)", "0.70", 1.0),
               (0.01, 0.02): ("r 0.01-0.02", "#92c5de", 1.2),
               (0.02, 0.02): ("r = 0.02", fs.INK, 2.0),
               (0.03, 0.03): ("r = 0.03", "#f4a582", 1.0),
               (0.04, 0.04): ("r = 0.04", fs.C_1984, 1.0)}
SOURCES = {"1996-start": ("2009-2014-1996", 1996), "2009-start": ("2009-2014", 2009)}
# -----------------------------------------------------------------------------


def save(fig, name, caption):
    path = OUT / TOPICS.get(name, "2-observed-dunes") / f"{name}.png"
    fs.save(fig, path, close=True)
    fs.record_caption(path, caption + " Producer: HAT_stable_dune_figures.py.")
    print(path, flush=True)


def finish(ax, i=None):
    ax.grid(axis="y", color=fs.GRID_C, lw=0.5)
    fs.open_frame(ax)
    if i is not None:
        fs._letter_inside(ax, i)


def shade_stable(ax):
    for label, doms in STABLE.items():
        ax.axvspan(min(doms) - 0.5, max(doms) + 0.5, color=STABLE_C[label], alpha=0.15, lw=0, zorder=0)


# --- data ----------------------------------------------------------------------------
def survey_table():
    t = pd.read_csv(tv.product_dir("1996-start") / "survey-crests" / "survey_dune_crests.csv")
    return t[t.status == "ok"]


# Domain-median crest per survey, on cells with a dune in every survey and not at the window's landward edge
def domain_crests():
    t = survey_table()
    w = t.pivot_table(index=["gis", "cell"], columns="survey", values="crest_m_navd88").dropna()
    edge = pd.read_csv(tv.product_dir("1996-start") / "survey-crests" / "survey_dune_crests.csv") \
        .groupby(["gis", "cell"]).at_landward_edge.any()
    w = w[~edge.reindex(w.index).fillna(True)]
    return w.groupby(level="gis").median(), w.groupby(level="gis").size()


def shoreline_lrr():
    d = pd.read_csv(obsr.domain_csv(1996, 2009))
    return d.set_index(d.domain_number.astype(int)).median_lrr


def overwash_count():
    o = pd.read_csv(ow.OBSERVATIONS / "overwash_observations.csv")
    o = o[(o.Year >= 1996) & (o.Year <= 2009)]
    return o.groupby("domain").overwash.sum()


# The 1996 start's cross-shore frame applied to any 10 m array on its grid (no sentinels): ocean first
def to_frame(a10, shear, c0):
    a = X.orient_ocean_right(a10, X.OCEAN_LOC)[:, ::-1]
    out = np.full(a.shape, np.nan)
    n = a.shape[1]
    for i in range(a.shape[0]):
        k = int(shear[i])
        out[i, :n - k] = a[i, k:]
    return out[:, c0:]


_WINDOWS = {}


# Elevation (m NAVD88) of one start in the 1996 frame, kept only where its own survey covers the whole 10 m cell
def profile_array(start, gis):
    if not _WINDOWS:
        _WINDOWS.update(json.loads(tv.picks_dir("1996-start").joinpath(
            f"HAT_dune_search_windows_{tv.topo_dirs('1996-start')[2]}.json").read_text(encoding="utf-8")))
    shear, c0, i0, i1 = S.frame("1996-start", gis, _WINDOWS)
    prod, year = SOURCES[start]
    root = product(prod).gapfill_1m
    with rasterio.open(root / f"clip_domain_{gis}_filled.tif") as s:
        z = S.clean(s.read(1), s.nodata)
    with rasterio.open(root / f"clip_domain_{gis}_survey.tif") as s:
        own = np.where(s.read(1) == year, 1.0, np.nan)
    z10 = S.downsample(z) * S.downsample(own)
    return to_frame(z10, shear, c0), (i0, i1)


# --- figures -------------------------------------------------------------------------
def fig_where():
    crest, n = domain_crests()
    chg = (crest[2009] - crest[1996]).where(n >= MIN_CELLS)
    lrr, ovw = shoreline_lrr(), overwash_count()
    gis = np.arange(1, 91)
    fig, (a, b, c) = plt.subplots(3, 1, figsize=fs.figsize("double", height=5.4), sharex=True,
                                  layout="constrained", gridspec_kw=dict(height_ratios=[1, 1, 0.55]))
    for ax in (a, b, c):
        shade_stable(ax)
    lr = lrr.reindex(gis)
    a.axhspan(-SHORE_BAND, SHORE_BAND, color="0.9", zorder=0)
    a.bar(gis, lr, 0.8, color=np.where(lr.abs() <= SHORE_BAND, fs.INK_MUTED, "0.78"))
    a.axhline(0, color=fs.INK, lw=0.6)
    a.set(ylabel="Shoreline change,\n1996–2009 (m/yr)")
    a.text(0.99, 0.95, f"grey band: within ±{SHORE_BAND} m/yr (stable)", transform=a.transAxes, ha="right",
           va="top", fontsize=6.5, color=fs.INK_MUTED)
    finish(a, 0)
    ch = chg.reindex(gis)
    b.bar(gis, ch, 0.8, color=np.where(ch.abs() <= 1.0, fs.INK_MUTED, "0.78"))
    b.axhline(0, color=fs.INK, lw=0.6)
    b.set(ylabel="Dune crest change,\n1996–2009 lidar (m)")
    b.text(0.99, 0.95, f"blank: fewer than {MIN_CELLS} cells with a dune in all nine surveys", transform=b.transAxes,
           ha="right", va="top", fontsize=6.5, color=fs.INK_MUTED)
    finish(b, 1)
    ov = ovw.reindex(gis).fillna(0)
    c.bar(gis, ov, 0.8, color=fs.C_1984)
    c.set(ylabel="Images with\noverwash", xlabel="Domain (GIS ID, south to north)", ylim=(0, max(3, ov.max() + 0.5)),
          xlim=(0.4, 90.6))
    c.set_xticks(range(0, 91, 10))
    finish(c, 2)
    for label, doms in STABLE.items():
        a.text(np.mean(doms), a.get_ylim()[1] * 0.97, label, ha="center", va="top", fontsize=7,
               color=STABLE_C[label], fontweight="bold")
    save(fig, "stable_dune_fig1_where",
         "Where the dunes were stable, 1996-2009. (a) Shoreline change rate per domain (CoastSat linear "
         "regression, 1996-2009); dark bars lie within +/-0.5 m/yr. (b) Change in the domain-median dune crest "
         "between the 1996 and 2009 lidar surveys, on the 10 m cells with a dune in all nine surveys. (c) Number of "
         "imagery dates showing overwash. Shaded: GIS 20-22 (Buxton-Avon / Avon) and GIS 59-62 (Mid-island / "
         "Wimble Shoals), where the shoreline held, no overwash was seen and the dune changed gently; these "
         "calibrate the growth rate.")


def profile_panel(ax, gis, i, note=None):
    xs = None
    for start, col, year in (("1996-start", C96, 1996), ("2009-start", C09, 2009)):
        z, (i0, i1) = profile_array(start, gis)
        med = np.nanmedian(z, axis=0)
        lo, hi = np.nanpercentile(z, 25, axis=0), np.nanpercentile(z, 75, axis=0)
        if xs is None:                                       # x = 0 at the 1996 shoreline (MHW crossing)
            j0 = int(np.argmax(med > X.MHW_M))
            xs = (np.arange(z.shape[1]) - j0) * fs.CELL_M
            win = ((i0 - j0) * fs.CELL_M, (i1 - j0) * fs.CELL_M)
        ax.fill_between(xs, lo, hi, color=col, alpha=0.18, lw=0)
        ax.plot(xs, med, color=col, lw=1.6, label=f"Observed: {year} lidar DEM")
    ax.axvspan(*win, color="0.93", zorder=0)
    ax.axhline(X.MHW_M, color=fs.INK_MUTED, lw=0.6, ls=":")
    ax.text(win[0] + 3, ax.get_ylim()[0] + 0.2 if False else -0.6, "dune window", fontsize=6, color=fs.INK_MUTED)
    ax.set(xlim=(-40, win[1] + 60), ylim=(-1, 9.5), xlabel="Distance landward of the 1996 shoreline (m)")
    title = f"GIS {gis}" + (f", {note}" if note else "")
    ax.text(0.97, 0.95, title, transform=ax.transAxes, ha="right", va="top", fontsize=7, fontweight="bold")
    finish(ax, i)


def fig_profiles(domains, name, caption, notes=None):
    fig, axes = plt.subplots(1, len(domains), figsize=fs.figsize("double", height=2.9), layout="constrained",
                             sharey=True)
    for i, (ax, g) in enumerate(zip(axes, domains)):
        profile_panel(ax, g, i, (notes or {}).get(g))
    axes[0].set_ylabel("Elevation (m NAVD88)")
    fig.legend(*axes[0].get_legend_handles_labels(), frameon=False, fontsize=6.5, loc="outside lower center", ncol=3)
    save(fig, name, caption)


def fig_surveys():
    crest, n = domain_crests()
    fig, axes = plt.subplots(1, 2, figsize=fs.figsize("double", height=2.9), layout="constrained", sharey=True)
    years = [int(y) for y in crest.columns]
    for i, (ax, (label, doms)) in enumerate(zip(axes, STABLE.items())):
        for g in doms:
            ax.plot(years, crest.loc[g].to_numpy(), color=STABLE_C[label], lw=0.8, alpha=0.5, marker="o", ms=2.5)
            ax.text(years[-1] + 0.3, crest.loc[g, 2009], f"{g}", fontsize=6, va="center", color=STABLE_C[label])
        ax.plot(years, crest.loc[list(doms)].median().to_numpy(), color=STABLE_C[label], lw=2, marker="o", ms=4,
                label=f"Observed: {label}, median")
        ax.axvline(2003.7, color=fs.INK_MUTED, lw=0.6, ls="--")
        ax.text(2003.8, 9.25, "Isabel", fontsize=6.5, color=fs.INK_MUTED)
        ax.axhline(5.58, color=fs.INK_MUTED, lw=0.6, ls="-.")
        ax.text(1996, 5.62, "Isabel run-up, 5.58 m", fontsize=6, color=fs.INK_MUTED)
        ax.set(xlabel="Survey year", xlim=(1995.5, 2010.2), ylim=(4, 9.5))
        ax.legend(frameon=False, fontsize=6.5, loc="lower right")
        finish(ax, i)
    axes[0].set_ylabel("Domain-median dune crest (m NAVD88)")
    save(fig, "stable_dune_fig4_surveys",
         "Dune crests at the stable stretches through all nine lidar surveys, 1996-2009 (median of the 10 m cells "
         "with a dune in every survey; thin lines: single domains, labelled). The crests sit above Isabel's run-up "
         "and change gradually, with no drop after Isabel (2004 survey), so the change between surveys is growth "
         "rather than storm damage: the process r controls. 1996-2005 surveys are 3-9 m grids, 2009 the 1 m DEM.")


def model_group_series(variant):
    mod = V.model_crests(R.START, variant)
    first = V.PERIODS[R.START][1]
    surv = V.surveys(R.START)
    out = {}
    for label, doms in STABLE.items():
        s = surv[surv.gis.isin(doms)]
        cells = s.dropna(subset=["crest"]).groupby(["gis", "cell"]).size()
        cells = cells[cells == s.survey.nunique()].index           # cells with a crest in every survey
        idx = np.array([(g - 1) * R.CELLS + c for g, c in cells])
        yrs = first + np.arange(mod.shape[0])
        out[label] = (yrs, np.median(mod[:, idx], axis=1), cells)
    return out


def fig_model():
    surv = V.surveys(R.START)
    fig, axes = plt.subplots(1, 2, figsize=fs.figsize("double", height=3.1), layout="constrained", sharey=True)
    series = {r: model_group_series(R.G.name(R.CEILING, r, groin=True)) for r in MODEL_RATES}
    for i, (ax, label) in enumerate(zip(axes, STABLE)):
        for r, (lab, col, lw) in MODEL_RATES.items():
            yrs, med, _ = series[r][label]
            ax.plot(yrs, med, color=col, lw=lw, label=f"Model, {lab}")
        cells = series[(0.0, 0.0)][label][2]
        s = surv.set_index(["gis", "cell"]).loc[cells].reset_index()
        pts = s.groupby("survey").agg(state=("model_state", "first"), crest=("crest", "median"))
        ax.plot([int(x[:4]) for x in pts.state], pts.crest, ls="none", marker="o", ms=5, mfc="white", mec=fs.INK,
                mew=1.3, label="Observed: lidar surveys", zorder=6)
        ax.text(0.97, 0.05, label, transform=ax.transAxes, ha="right", fontsize=7, fontweight="bold",
                color=STABLE_C[label])
        yrs, med, _ = series[(0.0, 0.0)][label]
        ax.annotate("year 1: road sand\nbulldozed onto dunes", (yrs[1], med[1]), xytext=(yrs[1] + 0.6, med[1] + 0.9),
                    fontsize=6, color=fs.INK_MUTED, arrowprops=dict(arrowstyle="-", color=fs.INK_MUTED, lw=0.6))
        ax.set(xlabel="Year (1 January)", xlim=(1996.5, 2010.5), ylim=(5, 9))
        finish(ax, i)
    axes[0].set_ylabel("Median dune crest (m NAVD88)")
    fig.legend(*axes[0].get_legend_handles_labels(), frameon=False, fontsize=6.5, loc="outside lower center", ncol=5)
    save(fig, "stable_dune_fig5_model",
         "Model against the surveys at the stable stretches, 1996 start, Dmax 10.5 m above the berm, groin on. "
         "Median crest of the cells with a dune in every survey, for four growth rates; circles: the lidar on the "
         "same cells, at the model state each survey's flight dates point to. With no growth the model dunes stay "
         "flat while the real ones rose; 0.02 follows the surveys on the corrected 1996 start (ALACE -0.24 m); "
         "faster rates outgrow them. The first-year rise at every rate, including no growth, is the roadway manager bulldozing "
         "3,400-5,100 m3 of sand above its road elevation onto the dunes once, in the road-managed domains "
         "(GIS 20 and 59-62; 21-22 are a community zone without it).")


# --- same dune in both DEMs ------------------------------------------------------------
# One year's 1 m DEM (m NAVD88) on the clip grid, NaN where that year's own survey does not cover the pixel
def dem_1m(start, gis):
    prod, year = SOURCES[start]
    root = product(prod).gapfill_1m
    with rasterio.open(root / f"clip_domain_{gis}_filled.tif") as s:
        z, t = S.clean(s.read(1), s.nodata), s.transform
    with rasterio.open(root / f"clip_domain_{gis}_survey.tif") as s:
        z[s.read(1) != year] = np.nan
    return z, t


# The crest the survey-crest code picks in each 10 m row (S.crests), with its position: the frame column, or -1
def crest_picks(start, gis):
    if not _WINDOWS:
        profile_array(start, gis)
    shear, c0, i0, i1 = S.frame("1996-start", gis, _WINDOWS)
    z1, _ = dem_1m(start, gis)
    z = S.profiles(S.downsample(z1), shear, c0)
    height, _, status = S.crests(z, i0, i1)
    above = z > X.BEACH_START_THR_M
    sb = np.where(above.any(axis=1), above.argmax(axis=1), -1)
    pos = np.full(z.shape[0], -1)
    for i in np.flatnonzero(status == "ok"):
        a = max(i0, int(sb[i])) if (X.CLIP_WINDOW_TO_BEACH and sb[i] >= 0) else i0
        win = z[i, a:min(i1, z.shape[1])]
        pos[i] = a + int(np.argmax(np.where(win > X.SENTINEL_WATER_M + 1e-9, win, -np.inf)))
    return pos, height + X.MHW_M, (shear, c0, i0, i1)


# A frame cell (row i, column j) back to UTM metres at its 10 m centre: undoes the flips, shear and water trim
def frame_to_xy(i, j, shear, c0, t, n_rows=50, n_cols=200):
    r, c = n_rows - 1 - i, n_cols - 1 - (j + c0 + shear[i])
    return t.c + 10 * c + 5, t.f - (10 * r + 5)


def fig_dem_maps(gis):
    z96, t = dem_1m("1996-start", gis)
    z09, _ = dem_1m("2009-start", gis)
    p96, h96, (shear, c0, i0, i1) = crest_picks("1996-start", gis)
    p09, h09, _ = crest_picks("2009-start", gis)
    rows = np.arange(len(p96))
    # Crop cross-shore to the beach and dune: 40 m seaward of the window to 80 m landward of it
    xs_lo = [frame_to_xy(i, i1 - 1, shear, c0, t)[0] for i in rows]
    xs_hi = [frame_to_xy(i, i0, shear, c0, t)[0] for i in rows]
    x0, x1 = min(xs_lo) - 80, max(xs_hi) + 40
    c_lo, c_hi = int(x0 - t.c), int(x1 - t.c)
    ext = (t.c + c_lo, t.c + c_hi, t.f - z96.shape[0], t.f)
    cmap, _, _ = fs.elevation_cmap()
    norm = matplotlib.colors.Normalize(0.0, 8.0)
    ls = matplotlib.colors.LightSource(azdeg=315, altdeg=45)
    diff = z09 - z96
    fig, axes = plt.subplots(1, 3, figsize=fs.figsize("double", height=5.2), layout="constrained",
                             sharex=True, sharey=True)
    for k, (ax, z, title) in enumerate(zip(axes, (z96, z09, diff),
                                           ("Observed: 1996 lidar (ALACE)", "Observed: 2009 lidar (USACE)",
                                            "2009 minus 1996"))):
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
        # The model's dune window: the band of 10 m cells each row's crest is searched in
        for j, side in ((i0, 1), (i1 - 1, -1)):
            xy = np.array([frame_to_xy(i, j, shear, c0, t) for i in rows])
            ax.step(xy[:, 0] + side * 5, xy[:, 1], where="mid", color=fs.INK, lw=0.7, ls="--")
        for pos, col, year in ((p96, C96, 1996), (p09, C09, 2009)):
            if (k == 0 and year == 2009) or (k == 1 and year == 1996):
                continue
            ok = pos >= 0
            xy = np.array([frame_to_xy(i, pos[i], shear, c0, t) for i in rows[ok]])
            ax.plot(xy[:, 0], xy[:, 1], ls="none", marker="o", ms=2.6, mfc=col, mec="white", mew=0.4,
                    label=f"Crest picked from the {year} lidar")
        ax.set_title(title, fontsize=7)
        ax.set_aspect("equal")
        ax.set(xlim=(ext[0], ext[1]), ylim=(ext[2], ext[3]), xticks=[], yticks=[])
        fs._letter_inside(ax, k)
    # Scale bar and north arrow on the first panel
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
    labels.append("Model dune window")
    fig.legend(handles, labels, frameon=False, fontsize=6.5, loc="outside lower center", ncol=3)
    both = (p96 >= 0) & (p09 >= 0)
    shift = (p09[both] - p96[both]) * fs.CELL_M
    label = next(lab for lab, doms in STABLE.items() if gis in doms)
    fig.suptitle(f"GIS {gis} ({label})", fontsize=8, fontweight="bold", x=0.01, ha="left")
    save(fig, f"stable_dune_fig7_dem_gis{gis}",
         f"The same dune in both surveys, GIS {gis}. (a) The 1996 lidar (ALACE, lowered 0.24 m) and (b) the 2009 "
         f"lidar (USACE) at 1 m, north up, true scale, over a hillshade; grey where that year's own survey has no "
         f"data (the 1996 start is the 2009 DEM there, so it is not shown). (c) 2009 minus 1996: blue gain, red "
         f"loss. Dashed: the model's dune window, the band of 10 m cells each row's crest is searched in. Dots: the "
         f"crest picked in each 10 m alongshore row from that year's DEM, averaged to 10 m as the model sees it. "
         f"The picks follow one ridge in both years: {np.mean(np.abs(shift) <= fs.CELL_M):.0%} of the "
         f"{both.sum()} rows picked in both years are within one 10 m cell of each other."
         + (" The white blocks landward of the dune in (a) are houses: the 1996 ALACE grid keeps buildings that "
            "the 2009 bare-earth DEM removes, hence the red patches in (c); they sit outside the dune window."
            if gis in STABLE["GIS 20-22"] else ""))


def fig_alongshore():
    t = survey_table()
    fig, axes = plt.subplots(2, 1, figsize=fs.figsize("double", height=4.2), layout="constrained", sharey=True)
    for i, (ax, (label, doms)) in enumerate(zip(axes, STABLE.items())):
        s = t[t.gis.isin(doms)].copy()
        s["x"] = (s.gis - min(doms)) * 500 + s.cell * fs.CELL_M + 5
        w = s.pivot_table(index="x", columns="survey", values="crest_m_navd88")
        for year in w.columns:
            if year not in (1996, 2009):
                ax.plot(w.index, w[year], color="0.8", lw=0.6, zorder=1)
        # 1996 and 2009 from the 1 m DEMs the maps show (the table's 1996 is NOAA's coarse grid, a bias check)
        for start, year, col in (("1996-start", 1996, C96), ("2009-start", 2009, C09)):
            h = np.concatenate([crest_picks(start, g)[1] for g in doms])
            w[year] = pd.Series(h, index=(np.arange(len(h)) * fs.CELL_M + 5)).reindex(w.index).to_numpy()
            ax.plot(w.index, w[year], color=col, lw=1.4, label=f"Observed: {year} lidar DEM (1 m)", zorder=3)
        both = w[[1996, 2009]].dropna()
        r = np.corrcoef(both[1996], both[2009])[0, 1]
        ax.plot([], [], color="0.8", lw=0.6, label="Observed: the 1997-2005 surveys (3-9 m grids)")
        for g in doms[1:]:
            ax.axvline((g - min(doms)) * 500, color=fs.INK_MUTED, lw=0.5, ls=":")
        for g in doms:
            ax.text((g - min(doms)) * 500 + 250, 9.3, f"GIS {g}", ha="center", fontsize=6.5, color=fs.INK_MUTED)
        ax.text(0.99, 0.04, f"{label}: 1996 vs 2009 alongshore r = {r:.2f}, median rise "
                f"{np.median(both[2009] - both[1996]):+.2f} m", transform=ax.transAxes, ha="right", fontsize=6.5,
                color=STABLE_C[label], fontweight="bold")
        ax.text(0.01, 0.04, label, transform=ax.transAxes, fontsize=7, fontweight="bold", color=STABLE_C[label])
        ax.set(xlim=(0, len(doms) * 500), ylim=(3, 9.8), ylabel="Dune crest (m NAVD88)")
        finish(ax, i)
    axes[-1].set_xlabel("Alongshore distance, south to north (m)")
    fig.legend(*axes[0].get_legend_handles_labels(), frameon=False, fontsize=6.5, loc="outside lower center", ncol=3)
    save(fig, "stable_dune_fig8_alongshore",
         "Dune crest along the stable stretches, 1996-2009, one value per 10 m model row (the highest 10 m cell in "
         "the model's dune window). Red and amber: the 1 m 1996 and 2009 DEMs the model's dunes and its 2009 target "
         "come from; grey: the 1997-2005 surveys, coarser grids. The 1996 and 2009 crests rise and dip at the same places "
         "alongshore, so the window picks the same dune in both surveys, and the 2009 line sits a little above the "
         "1996 one almost everywhere: the dune grew in place. Gaps: rows a survey did not cover. The two single-row 1996 spikes in GIS 21 and 22 "
         "(above 9 m) appear in no other survey.")


def fig_crest_shift():
    rows = []
    for label, doms in STABLE.items():
        for g in doms:
            p96, _, _ = crest_picks("1996-start", g)
            p09, _, _ = crest_picks("2009-start", g)
            ok = (p96 >= 0) & (p09 >= 0)
            rows += [dict(group=label, shift=s) for s in (p09[ok] - p96[ok]) * fs.CELL_M]
    d = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=fs.figsize("single", height=2.7), layout="constrained")
    bins = np.arange(-65, 70, 10)
    ax.hist([d[d.group == lab]["shift"] for lab in STABLE], bins=bins, stacked=True,
            color=[STABLE_C[lab] for lab in STABLE], label=[f"Observed: {lab}" for lab in STABLE], rwidth=0.85)
    within = np.mean(np.abs(d["shift"]) <= fs.CELL_M)
    ax.axvline(0, color=fs.INK, lw=0.6, ls="--")
    ax.text(0.98, 0.95, f"{within:.0%} of {len(d)} rows within\none 10 m cell", transform=ax.transAxes,
            ha="right", va="top", fontsize=6.5)
    ax.set(xlabel="Crest position, 2009 minus 1996 (m, + landward)", ylabel="10 m rows")
    ax.legend(frameon=False, fontsize=6.5, loc="upper left")
    finish(ax)
    save(fig, "stable_dune_fig9_crest_shift",
         "How far the picked dune crest moved cross-shore between the 1996 and 2009 lidar, per 10 m row at the "
         "seven stable domains (rows picked in both years). Positions are 10 m model cells, so 0 means the same "
         f"cell. {within:.0%} of rows are within one cell: the crest picked in 2009 is the 1996 dune, not a "
         "different ridge.")


def main():
    fs.apply_style()
    OUT.mkdir(parents=True, exist_ok=True)
    if len(sys.argv) > 1 and sys.argv[1] == "dems":           # only the same-dune figures (no model runs needed)
        for g in sorted(g for doms in STABLE.values() for g in doms):
            fig_dem_maps(g)
        fig_alongshore()
        fig_crest_shift()
        return
    fig_where()
    fig_profiles(PROFILE_DOMAINS, "stable_dune_fig2_profiles",
                 "Cross-shore profiles at the two stable stretches from the 1996 and 2009 lidar DEMs, in the same "
                 "frame (the 1996 start's, 10 m cells; line: median of the domain's 50 alongshore rows, band: "
                 "25th-75th percentile; only cells covered by that year's own survey). Grey: the dune window the "
                 "model reads its dune from. Shoreline and dune stayed in place and the crest rose slightly: a "
                 "dune that grew without storm loss or shoreline change.")
    fig_profiles(tuple(EXCLUDED), "stable_dune_fig3_excluded",
                 "As stable_dune_fig2 at two domains left out of the calibration. GIS 44: the beach built out "
                 "(+1.7 m/yr) and the dune grew 2 m, fed by the extra sand, which would overstate r. GIS 65: the "
                 "shoreline retreated (-2.2 m/yr); the crest held while the dune lost width, which r cannot "
                 "represent.", EXCLUDED)
    fig_surveys()
    fig_model()


if __name__ == "__main__":
    main()
