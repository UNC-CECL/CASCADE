"""
Why the stable dunes at GIS 59-62 and 20-22 can calibrate the dune growth rate r, one point per figure.

    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/dunes/HAT_stable_dune_figures.py

Figures for the stable-dune study (HAT_stable_dune_growth_rate.py, run it first),
written to its figures/ folder (data/.../1996-start/dune-growth-rate/ since 2026-10-21), one subfolder per
topic: 1-selection (fig1, fig3), 2-observed-dunes, 3-model-fit (the fit figures: stable-stretches/ for why
0.01-0.03, island-wide/ for how it does across the island, 2026-10-22). 2-observed-dunes is split by
domain: all-stable-domains/ (fig2, fig4, fig8, fig9) and one folder per stable domain, gis<N>/, with every
observed-dune figure for that domain alone (stable_dune_gis<N>_*). In reading order:
  stable_dune_fig1_where       shoreline change, dune change and overwash per domain, 1996-2009:
                               the two stretches meet all three conditions
  stable_dune_fig2_profiles    1996 and 2009 DEM profiles at the stable stretches: dune and beach in place
  stable_dune_fig3_excluded    the same at two excluded domains: a growing beach, a retreating shoreline
  stable_dune_fig4_surveys     crests through all nine 1996-2009 lidar surveys: no storm losses
  stable_dune_gis<N>_dem       the 1996 and 2009 1 m DEMs and their difference, north up, with the dune window
                               and each row's picked crest: the same dune in both (was fig7)
  stable_dune_gis<N>_profile, _surveys, _alongshore, _crest_shift
                               fig2, fig4, fig8 and fig9 for that domain alone
  stable_dune_fig8_alongshore  crest along the shore in every survey: 1996 and 2009 rise and dip together
  stable_dune_fig9_crest_shift how far the picked crest moved cross-shore, 1996 to 2009
    (the gis<N> sets, fig8 and fig9 alone, no model runs needed: HAT_stable_dune_figures.py dems)
  stable_dune_fit_trajectory   why 0.01-0.03: the adopted rate against every survey, cell spread shown
  stable_dune_fit_error_by_rate every rate tested, from tables/stable_dune_summary.csv
  stable_dune_fit_spread       a range keeps the median and doubles the spread of a fixed rate
  stable_dune_fit_island       the adopted rate island-wide: survey medians and low dunes through 1996-2009
  stable_dune_fit_island_alongshore  domain-median crest along the island, 2010 model vs 2009 lidar: (a) the
                               full-management run, (b) the natural run, groin on; each with the domain-mean
                               fit under the median one (was four figures: _natural, _mean, _natural_mean)
  stable_dune_fit_island_alongshore_1996_2009_mean.gif  the full-management alongshore fit animated with domain
                               means: one frame per model year, each survey on its scoring year
                               (HAT_stable_dune_figures.py gif)
    (the fit figures alone: HAT_stable_dune_figures.py fit; the alongshore figure alone: ... alongshore)
Profiles: the 1 m gap-filled DEMs of the 1996 and 2009 starts (same clip grid),
averaged to 10 m as the starts are and put in the 1996 start's cross-shore frame
by the survey-crest code. A 10 m cell is drawn only where all its core 1 m
pixels come from that year's own survey (1996 ALACE, 2009 USACE); the 1996 start
is the 2009 DEM wherever ALACE has no data, so unmasked cells would fake stability.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-26
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
sys.path.insert(0, str(next(_p for _p in _HERE.parents if (_p / "experiment_paths.py").exists())))
import experiment_paths  # noqa: E402,F401  (every experiments/<theme>/ folder onto sys.path)
# Hatteras/ holds the scripts; data/ and output/ live in the workspace on D: (hatteras_paths)
import sys as _hat_sys  # noqa: E402
_HAT_DIR = next(p for p in _HERE.parents  if (p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE as _REPO  # noqa: E402
sys.path.insert(0, str(_HAT_DIR / "input_prep" / "1-barrier3d-domains" / "2-survey-crests"))
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
# One subfolder per point the figures make, in reading order (2026-10-21); 3-model-fit split 2026-10-22
TOPICS = {"stable_dune_fig1_where": "1-selection", "stable_dune_fig3_excluded": "1-selection",
          "stable_dune_fit_trajectory": "3-model-fit/stable-stretches",
          "stable_dune_fit_error_by_rate": "3-model-fit/stable-stretches",
          "stable_dune_fit_spread": "3-model-fit/stable-stretches",
          "stable_dune_fit_island": "3-model-fit/island-wide",
          "stable_dune_fit_island_alongshore": "3-model-fit/island-wide"}
# Everything else is 2-observed-dunes, split by the domains a figure covers: one folder per stable domain
# (stable_dune_gis<N>_*) and all-stable-domains for the figures that compare them (2026-10-21)
OBSERVED = "2-observed-dunes"
C96, C09 = fs.C_1984, "#c8880f"            # house period colours: 1996 red, 2009 amber
STABLE = R.GROUPS                           # {"GIS 59-62": (...), "GIS 20-22": (...)}
STABLE_C = R.GROUP_COLOURS
PROFILE_DOMAINS = (61, 21)                  # one per stable stretch
EXCLUDED = {44: "beach building out", 65: "shoreline retreating"}
SHORE_BAND = 0.5                            # m/yr: the stable-shoreline condition
MIN_CELLS = 25                              # domains with fewer fully surveyed cells are left blank in fig1
SOURCES = {"1996-start": ("2009-2014-1996", 1996), "2009-start": ("2009-2014", 2009)}
# The current runner's 1996 full-management run (HAT_run_matrix matrix, 2026-10-21: spread groin 0.10, ends re-solved),
# for the island-wide alongshore figure; the study's own runs carry the 2026-10-20 groin and ends
MATRIX_RUN = R.U.PROJECT_ROOT / "output" / "raw_runs" / "matrix" / "1996_2009" / "edgeBE" \
    / "HAT_1996_2009_edgeBE_offsetmetres_road_bdm_groin"
# The same matrix's natural run (no roadway or beach/dune management, groin on), for comparison
MATRIX_RUN_NATURAL = MATRIX_RUN.with_name("HAT_1996_2009_edgeBE_offsetmetres_noroad_nobdm_groin")
# -----------------------------------------------------------------------------


def folder(name):
    if name in TOPICS:
        return TOPICS[name]
    if name.startswith("stable_dune_gis"):
        return f"{OBSERVED}/gis{name.split('_')[2][3:]}"
    return f"{OBSERVED}/all-stable-domains"


def save(fig, name, caption):
    path = OUT / folder(name) / f"{name}.png"
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


# Crest per cell and survey, on cells with a dune in every survey and not at the window's landward edge
def cell_crests():
    t = survey_table()
    w = t.pivot_table(index=["gis", "cell"], columns="survey", values="crest_m_navd88").dropna()
    edge = pd.read_csv(tv.product_dir("1996-start") / "survey-crests" / "survey_dune_crests.csv") \
        .groupby(["gis", "cell"]).at_landward_edge.any()
    return w[~edge.reindex(w.index).fillna(True)]


def domain_crests():
    w = cell_crests()
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
         "between the 1996 and 2009 lidar surveys, on the 10 m cells with a dune in all nine surveys; blank where a "
         "domain has fewer than 25 such cells. (c) Number of "
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


# Model crest per saved state, GIS 1-90 x 50 cells, m NAVD88. The run's own buffer count, not today's config:
# these runs have 15 buffers a side, the site config 30 since 2026-10-21, so DOM.gis_to_pad would read the wrong
# domains. Buffers are symmetric, so GIS 1 sits (domains - 90) / 2 in
_CRESTS = {}


def model_crests(variant):
    if variant not in _CRESTS:
        c = R.U.load(variant if isinstance(variant, Path) else R.G.run_dir(R.START, variant))
        pad = (len(c.barrier3d) - 90) // 2
        cc = [np.asarray(c.barrier3d[pad + g - 1].DuneDomain, float).max(axis=2) * 10 for g in range(1, 91)]
        _CRESTS[variant] = np.stack(cc, axis=1).reshape(cc[0].shape[0], -1) + V.BERM_NAVD88
    return _CRESTS[variant]


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
    save(fig, f"stable_dune_gis{gis}_dem",
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


# One stretch (or one domain) of crests alongshore; returns the 1996-2009 alongshore r and median rise
def alongshore_panel(ax, doms, label, i=None, name=None):
    t = survey_table()
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
    for g in list(doms)[1:]:
        ax.axvline((g - min(doms)) * 500, color=fs.INK_MUTED, lw=0.5, ls=":")
    for g in doms if name is None else ():        # one domain: its name is in the corner label
        ax.text((g - min(doms)) * 500 + 250, 9.3, f"GIS {g}", ha="center", fontsize=6.5, color=fs.INK_MUTED)
    ax.text(0.99, 0.04, f"{name or label}: 1996 vs 2009 alongshore r = {r:.2f}, median rise "
            f"{np.median(both[2009] - both[1996]):+.2f} m", transform=ax.transAxes, ha="right", fontsize=6.5,
            color=STABLE_C[label], fontweight="bold")
    if name is None:
        ax.text(0.01, 0.04, label, transform=ax.transAxes, fontsize=7, fontweight="bold", color=STABLE_C[label])
    ax.set(xlim=(0, len(doms) * 500), ylim=(3, 9.8), ylabel="Dune crest (m NAVD88)")
    finish(ax, i)
    return r, np.median(both[2009] - both[1996])


def fig_alongshore():
    fig, axes = plt.subplots(2, 1, figsize=fs.figsize("double", height=4.2), layout="constrained", sharey=True)
    for i, (ax, (label, doms)) in enumerate(zip(axes, STABLE.items())):
        alongshore_panel(ax, doms, label, i)
    axes[-1].set_xlabel("Alongshore distance, south to north (m)")
    fig.legend(*axes[0].get_legend_handles_labels(), frameon=False, fontsize=6.5, loc="outside lower center", ncol=3)
    save(fig, "stable_dune_fig8_alongshore",
         "Dune crest along the stable stretches, 1996-2009, one value per 10 m model row (the highest 10 m cell in "
         "the model's dune window). Red and amber: the 1 m 1996 and 2009 DEMs the model's dunes and its 2009 target "
         "come from; grey: the 1997-2005 surveys, coarser grids. The 1996 and 2009 crests rise and dip at the same places "
         "alongshore, so the window picks the same dune in both surveys, and the 2009 line sits a little above the "
         "1996 one almost everywhere: the dune grew in place. Gaps: rows a survey did not cover. The two single-row 1996 spikes in GIS 21 and 22 "
         "(above 9 m) appear in no other survey.")


# Cross-shore move of the picked crest, 2009 minus 1996 (m, + landward), on rows picked in both years
def crest_shift(g):
    p96, _, _ = crest_picks("1996-start", g)
    p09, _, _ = crest_picks("2009-start", g)
    ok = (p96 >= 0) & (p09 >= 0)
    return (p09[ok] - p96[ok]) * fs.CELL_M


def fig_crest_shift():
    rows = []
    for label, doms in STABLE.items():
        for g in doms:
            rows += [dict(group=label, shift=s) for s in crest_shift(g)]
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


# --- one stable domain: every observed-dune figure for it, in its own folder -------------
def stretch_of(gis):
    return next(lab for lab, doms in STABLE.items() if gis in doms)


def fig_domain_profile(gis):
    fig, ax = plt.subplots(figsize=fs.figsize("single", height=2.8), layout="constrained")
    profile_panel(ax, gis, None)
    ax.set_ylabel("Elevation (m NAVD88)")
    ax.legend(frameon=False, fontsize=6.5, loc="upper left")
    save(fig, f"stable_dune_gis{gis}_profile",
         f"Cross-shore profile at GIS {gis} (in the {stretch_of(gis)} stretch) from the 1996 and 2009 lidar DEMs, in the same "
         f"frame (the 1996 start's, 10 m cells; line: median of the domain's 50 alongshore rows, band: 25th-75th "
         f"percentile; only cells covered by that year's own survey). Grey: the dune window the model reads its "
         f"dune from. As stable_dune_fig2_profiles, for this domain alone.")


def fig_domain_surveys(gis):
    w = cell_crests().loc[gis]
    years = [int(y) for y in w.columns]
    lab = stretch_of(gis)
    fig, ax = plt.subplots(figsize=fs.figsize("single", height=2.8), layout="constrained")
    ax.fill_between(years, w.quantile(0.25).to_numpy(), w.quantile(0.75).to_numpy(), color=STABLE_C[lab],
                    alpha=0.2, lw=0, label="Observed: middle half of the cells")
    ax.plot(years, w.median().to_numpy(), color=STABLE_C[lab], lw=2, marker="o", ms=4,
            label=f"Observed: GIS {gis}, median")
    ax.axvline(2003.7, color=fs.INK_MUTED, lw=0.6, ls="--")
    ax.text(2003.8, 9.25, "Isabel", fontsize=6.5, color=fs.INK_MUTED)
    ax.axhline(5.58, color=fs.INK_MUTED, lw=0.6, ls="-.")
    ax.text(1996, 5.62, "Isabel run-up, 5.58 m", fontsize=6, color=fs.INK_MUTED)
    ax.set(xlabel="Survey year", ylabel="Dune crest (m NAVD88)", xlim=(1995.5, 2010.2), ylim=(4, 9.5))
    ax.legend(frameon=False, fontsize=6.5, loc="lower right")
    finish(ax)
    rise = w[2009].median() - w[1996].median()
    save(fig, f"stable_dune_gis{gis}_surveys",
         f"Dune crest at GIS {gis} (in the {lab} stretch) through all nine lidar surveys, 1996-2009: median and middle half of "
         f"the {len(w)} 10 m cells with a dune in every survey. The median rose {rise:+.2f} m from 1996 to 2009, "
         f"with no drop after Isabel (2004 survey). 1996-2005 surveys are 3-9 m grids, 2009 the 1 m DEM. As "
         f"stable_dune_fig4_surveys, for this domain alone.")


def fig_domain_alongshore(gis):
    lab = stretch_of(gis)
    fig, ax = plt.subplots(figsize=fs.figsize("double", height=2.5), layout="constrained")
    r, rise = alongshore_panel(ax, (gis,), lab, name=f"GIS {gis}")
    ax.set_xlabel("Alongshore distance, south to north (m)")
    fig.legend(*ax.get_legend_handles_labels(), frameon=False, fontsize=6.5, loc="outside lower center", ncol=3)
    save(fig, f"stable_dune_gis{gis}_alongshore",
         f"Dune crest along GIS {gis} (in the {lab} stretch), one value per 10 m model row (the highest 10 m cell in the model's "
         f"dune window). Red and amber: the 1 m 1996 and 2009 DEMs; grey: the 1997-2005 surveys, coarser grids. "
         f"1996 and 2009 correlate alongshore at r = {r:.2f} and the median row rose {rise:+.2f} m. Gaps: rows a "
         f"survey did not cover. As stable_dune_fig8_alongshore, for this domain alone." + spike_note(gis))


# A single 1996 row above 9 m, seen in no other survey, drags r down at GIS 21 and 22: report r without it
def spike_note(gis):
    a, b = crest_picks("1996-start", gis)[1], crest_picks("2009-start", gis)[1]
    ok = ~np.isnan(a) & ~np.isnan(b)
    spike = ok & (a >= 9)
    if not spike.any():
        return ""
    keep = ok & ~spike
    return (f" The single-row 1996 spike above 9 m appears in no other survey; without it r = "
            f"{np.corrcoef(a[keep], b[keep])[0, 1]:.2f}.")


def fig_domain_crest_shift(gis):
    d = crest_shift(gis)
    lab = stretch_of(gis)
    fig, ax = plt.subplots(figsize=fs.figsize("single", height=2.6), layout="constrained")
    ax.hist(d, bins=np.arange(-65, 70, 10), color=STABLE_C[lab], rwidth=0.85, label=f"Observed: GIS {gis}")
    within = np.mean(np.abs(d) <= fs.CELL_M)
    ax.axvline(0, color=fs.INK, lw=0.6, ls="--")
    ax.text(0.98, 0.95, f"{within:.0%} of {len(d)} rows within\none 10 m cell", transform=ax.transAxes,
            ha="right", va="top", fontsize=6.5)
    ax.set(xlabel="Crest position, 2009 minus 1996 (m, + landward)", ylabel="10 m rows")
    ax.legend(frameon=False, fontsize=6.5, loc="upper left")
    finish(ax)
    save(fig, f"stable_dune_gis{gis}_crest_shift",
         f"How far the picked dune crest moved cross-shore between the 1996 and 2009 lidar at GIS {gis} (in the {lab} stretch), "
         f"per 10 m row picked in both years. 0 means the same 10 m cell; {within:.0%} of rows are within one "
         f"cell. As stable_dune_fig9_crest_shift, for this domain alone.")


def fig_domain(gis):
    fig_dem_maps(gis)
    fig_domain_profile(gis)
    fig_domain_surveys(gis)
    fig_domain_alongshore(gis)
    fig_domain_crest_shift(gis)


# --- why 0.01-0.03 ---------------------------------------------------------------------
ADOPTED = (0.01, 0.03)
SPREAD_RATES = {(0.02, 0.02): "Model, r = 0.02 (fixed)", (0.01, 0.03): "Model, r 0.01-0.03 (adopted)",
                (0.005, 0.035): "Model, r 0.005-0.035"}


def rate_label(r):
    return f"r = {r[0]:g}" if r[0] == r[1] else f"r {r[0]:g}-{r[1]:g}"


# The stable cells with a crest in every survey, per stretch, as flat model indices; and the survey table on them
def stable_cells():
    surv = V.surveys(R.START)
    out = {}
    for label, doms in STABLE.items():
        s = surv[surv.gis.isin(doms)]
        n = s.dropna(subset=["crest"]).groupby(["gis", "cell"]).size()
        cells = n[n == s.survey.nunique()].index
        out[label] = (np.array([(g - 1) * R.CELLS + c for g, c in cells]),
                      s.set_index(["gis", "cell"]).loc[cells].reset_index())
    return out


def model_at(mod, state):
    return mod[int(state[:4]) - V.PERIODS[R.START][1]]


def fig_fit_trajectory():
    cells = stable_cells()
    runs = {r: model_crests(R.G.name(R.CEILING, r, groin=True)) for r in ((0.0, 0.0), ADOPTED, (0.04, 0.04))}
    fig, axes = plt.subplots(1, 2, figsize=fs.figsize("double", height=3.0), layout="constrained", sharey=True)
    for i, (ax, label) in enumerate(zip(axes, STABLE)):
        idx, s = cells[label]
        yrs = V.PERIODS[R.START][1] + np.arange(runs[ADOPTED].shape[0])
        m = runs[ADOPTED][:, idx]
        ax.fill_between(yrs, *np.percentile(m, [25, 75], axis=1), color=fs.INK, alpha=0.12, lw=0,
                        label="Model, r 0.01-0.03: middle half of cells")
        ax.plot(yrs, np.median(m, axis=1), color=fs.INK, lw=2.0, label="Model, r 0.01-0.03 (adopted): median")
        for r, ls in (((0.0, 0.0), ":"), ((0.04, 0.04), "--")):
            med = np.median(runs[r][:, idx], axis=1)
            ax.plot(yrs, med, color=fs.INK_MUTED, lw=0.9, ls=ls)
            ax.text(yrs[-1] + 0.15, med[-1], "no growth" if r[0] == 0 else "r = 0.04", fontsize=6, va="center",
                    color=fs.INK_MUTED)
        g = s.groupby("survey").agg(state=("model_state", "first"), med=("crest", "median"),
                                    q1=("crest", lambda x: x.quantile(0.25)), q3=("crest", lambda x: x.quantile(0.75)))
        x = np.array([int(v[:4]) for v in g.state]) + 0.08     # nudge off the model line so both show
        ax.errorbar(x, g.med, yerr=[g.med - g.q1, g.q3 - g.med], fmt="o", ms=5, mfc="white", mec=STABLE_C[label],
                    mew=1.4, ecolor=STABLE_C[label], elinewidth=1, capsize=0, zorder=6,
                    label="Observed: lidar survey, median and middle half")
        ax.text(0.97, 0.05, f"{label} ({len(idx)} cells)", transform=ax.transAxes, fontsize=7, fontweight="bold",
                color=STABLE_C[label], ha="right")
        ax.set(xlabel="Year (1 January)", xlim=(1996.6, 2011.3))
        finish(ax, i)
    axes[0].set_ylabel("Dune crest (m NAVD88)")
    fig.legend(*axes[0].get_legend_handles_labels(), frameon=False, fontsize=6.5, loc="outside lower center", ncol=3)
    save(fig, "stable_dune_fit_trajectory",
         "The adopted growth rate follows the stable dunes through 1996-2009. Model crests with r 0.01-0.03 (line: "
         "median; band: middle half of the cells) against every lidar survey on the same cells (circles: median; "
         "bars: middle half), at the model state each survey's flight dates point to. 1996 start corrected for the "
         "ALACE offset, Dmax 10.5 m above the berm, groin on. Dotted: no growth; dashed: r = 0.04. 1996-2005 "
         "surveys are 3-9 m grids, which read sharp crests somewhat low. The first-year rise in every run, no growth "
         "included, is the roadway manager bulldozing 3,400-5,100 m3 of sand above its road elevation onto the dunes "
         "once, in the road-managed domains (GIS 20 and 59-62; 21-22 are a community zone without it).")


# The adopted rate island-wide, from the 1996 run with r 0.01-0.03: the model crests and the surveys
def island_data(run=None):
    mod = model_crests(run or R.G.name(R.CEILING, ADOPTED, groin=True))
    surv = V.surveys(R.START).dropna(subset=["crest"])
    return mod, surv, V.PERIODS[R.START][1]


# Median and 10th-percentile crest through the surveys, on the cells each survey covers (as
# HAT_dune_growth_low_rates.trajectory scores it); the model lines on the cells of the 2009 lidar
def fig_fit_island():
    mod, surv, first = island_data(MATRIX_RUN)
    yrs = first + np.arange(mod.shape[0])
    fig, ax = plt.subplots(figsize=fs.figsize("single", height=3.0), layout="constrained")
    rows = []
    for year, g in surv.groupby("survey"):
        k = int(g.model_state.iloc[0][:4]) - first
        if k <= 0 or k >= mod.shape[0]:
            continue
        idx = ((g.gis - 1) * R.CELLS + g.cell).to_numpy()
        rows.append(dict(x=first + k, obs_med=g.crest.median(), obs_p10=g.crest.quantile(0.1),
                         mod_med=np.median(mod[k][idx]), mod_p10=np.percentile(mod[k][idx], 10)))
    t = pd.DataFrame(rows)
    on09 = ((surv[surv.survey == 2009].gis - 1) * R.CELLS + surv[surv.survey == 2009].cell).to_numpy()
    for stat, ls, name in (("med", "-", "median"), ("p10", "--", "10th percentile")):
        ax.plot(yrs, np.percentile(mod[:, on09], 50 if stat == "med" else 10, axis=1), color=fs.INK, lw=1.6, ls=ls,
                label=f"Model, r 0.01-0.03: {name}")
        ax.plot(t.x + 0.08, t[f"obs_{stat}"], "o" if stat == "med" else "s", ms=5, mfc="white", mec=fs.C_1984,
                mew=1.4, ls="none", zorder=6, label=f"Observed: lidar survey {name}")
    end = t.iloc[-1]
    ax.text(0.97, 0.72, f"2009 median: model {end.mod_med:.2f}, lidar {end.obs_med:.2f} m", transform=ax.transAxes,
            ha="right", fontsize=6.5, color=fs.INK)
    ax.set(xlabel="Year (1 January)", ylabel="Dune crest (m NAVD88)", xlim=(1996.6, 2010.6))
    finish(ax)
    fig.legend(*ax.get_legend_handles_labels(), frameon=False, fontsize=6.5, loc="outside lower center", ncol=2)
    save(fig, "stable_dune_fit_island",
         "The adopted growth rate island-wide, 1996-2009. Median (circles, solid) and 10th-percentile (squares, "
         "dashed) dune crest over all surveyed cells, model with r 0.01-0.03 against every lidar survey on the cells "
         "it covers, at the model state its flight dates point to. Model: the current runner's 1996 full-management "
         "run (matrix, 2026-10-21): 1996 start corrected for the ALACE offset, Dmax 10.5 m above the berm, spread "
         "groin 0.10, edgeBE ends re-solved, storms v3_240, option A waves. The model's 10th percentile sits high "
         "through 2006: low dunes eroded 1996-2004 and stayed low in the surveys, while the model holds them up "
         "until 2006 and only then drops to the surveyed low by 2009.")


# Domain-median crest along the island, 1 January 2010 model vs the 2009 lidar; villages, shoals, groin and piers
# from the site config, as the other alongshore figures draw them. One panel per run (full management, natural);
# the domain-mean fit is printed under the median one (2026-10-22 check: the median choice barely moves it), so
# the four separate _natural/_mean figures are one figure since 2026-10-26
SCENARIOS = {"full management": (MATRIX_RUN, "Full management, groin on"),
             "natural, groin on": (MATRIX_RUN_NATURAL, "Natural (no roadway or beach/dune management), groin on")}


def island_fit_stats(obs, m):
    d = (m.reindex(obs.index) - obs).dropna()
    return np.corrcoef(obs[d.index], m[d.index])[0, 1], np.sqrt((d ** 2).mean()), d.mean(), len(d)


def alongshore_fit_panel(ax, run, title, i):
    from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS as A
    from cascade_pipeline.annotations import AnnotationConfig
    mod, surv, first = island_data(run)
    s09 = surv[surv.survey == 2009]
    s09 = s09.assign(m=mod[2010 - first][((s09.gis - 1) * R.CELLS + s09.cell).to_numpy()])
    stats = {stat: island_fit_stats(s09.groupby("gis").crest.agg(stat), s09.groupby("gis").m.agg(stat))
             for stat in ("median", "mean")}
    obs, m10 = s09.groupby("gis").crest.median(), s09.groupby("gis").m.median()
    shoal_c = AnnotationConfig.color_shoal
    for name, (lo, hi) in A.shoal_zones.items():
        ax.axvspan(lo - 0.5, hi + 0.5, facecolor=shoal_c, alpha=0.10, hatch="///", edgecolor=shoal_c, lw=0, zorder=0)
        ax.text((lo + hi) / 2, 0.03, name, transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                fontsize=6.5, color="#7A5800", style="italic", zorder=7, path_effects=fs._halo(2.0))
    ax.plot(m10.index, m10, color=fs.INK, lw=1.4, zorder=4, label="Model, 1 January 2010")
    ax.plot(obs.index, obs, color=C09, lw=0.9, marker="o", ms=2.8, zorder=5, label="Observed, 2009 lidar")
    ax.set(ylabel="Domain-median dune\ncrest (m NAVD88)", xlim=(0.4, 90.6), ylim=(0, 10.6))
    ax.set_xticks(range(0, 91, 10))
    # Above the axes at the right, clear of the shoal and structure labels: the plotted (median) fit, the
    # domain-mean fit under it in grey
    for k, (stat, colour) in enumerate((("median", fs.INK), ("mean", fs.INK_MUTED))):
        r, rmse, bias, _ = stats[stat]
        ax.text(1.0, 1.02 + 0.06 * (1 - k), f"domain {stat}:  r = {r:.2f}   RMSE = {rmse:.2f} m   "
                f"bias = {bias:+.2f} m", transform=ax.transAxes, ha="right", va="bottom",
                fontsize=6.5, color=colour)
    ax.set_title(f"{fs._letter(i)}  {title}", loc="left", fontweight="bold")
    finish(ax)
    fs.town_bands(ax, where="top", strip=0.07)
    return stats


def fig_fit_island_alongshore():
    from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS as A
    fig, axes = plt.subplots(2, 1, figsize=fs.figsize("double", height=5.6), layout="constrained", sharex=True)
    stats = {sc: alongshore_fit_panel(ax, run, title, i)
             for i, (ax, (sc, (run, title))) in enumerate(zip(axes, SCENARIOS.items()))}
    axes[-1].set_xlabel("Domain (GIS ID)")
    for x, ha, text in ((0.0, "left", A.low_end_label), (1.0, "right", A.high_end_label)):
        axes[-1].annotate(text, (x, 0), xycoords="axes fraction", xytext=(0, -24), textcoords="offset points",
                          ha=ha, va="top", fontsize=6.5, color=fs.INK_MUTED)
    fig.legend(*axes[0].get_legend_handles_labels(), frameon=False, fontsize=6.5, loc="outside lower center", ncol=2)
    for ax in axes:
        fs.structures(ax)
    fit = "; ".join(f"{letter} {sc}: r = {stats[sc]['median'][0]:.2f}, RMSE {stats[sc]['median'][1]:.2f} m, "
                    f"bias {stats[sc]['median'][2]:+.2f} m (domain means: r = {stats[sc]['mean'][0]:.2f}, "
                    f"RMSE {stats[sc]['mean'][1]:.2f} m, bias {stats[sc]['mean'][2]:+.2f} m)"
                    for letter, sc in zip(("(a)", "(b)"), SCENARIOS))
    save(fig, "stable_dune_fit_island_alongshore",
         f"Domain-median dune crest along Hatteras Island, the 1 January 2010 model with r 0.01-0.03 against the "
         f"2009 lidar, over the 10 m dune cells the lidar covers, across {stats['full management']['median'][3]} "
         f"domains. (a) the full-management run; (b) the natural run (no roadway or beach/dune management), groin "
         f"on. {fit}. The grey line of numbers in each panel is the same fit with domain means, which changes "
         f"little. Top strip: villages; hatched: Avon and Wimble Shoals; solid line: Buxton groin; dotted: piers. "
         f"Model: the current runner's 1996 matrix runs (2026-10-21): 1996 start corrected for the ALACE offset, "
         f"Dmax 10.5 m above the berm, spread groin 0.10, edgeBE ends re-solved, storms v3_240, option A waves.")


# fig_fit_island_alongshore through time, one frame per model state (1 January 1997-2010), kept bare: model line,
# survey, date. The model line is on the 2009 lidar's cells in every frame so it moves only by growth; a survey
# appears on the frame its scoring pairs it with (1 January the year after) and stays, faint, until the next.
# 1999 and 2000 both pair with 1 January 2000: the later one is drawn. Stats on that survey's own cells.
# Domain means need the lidar spikes out of the survey crests: HAT_survey_dune_crests.py blanks them (status
# "spike", 2026-10-26), so V.surveys no longer carries them.
GIF_FRAME_S, GIF_HOLD_S, GIF_DPI = 0.8, 3.0, 130


def gif_fit_island_alongshore(stat="median"):
    from PIL import Image
    mod, surv, first = island_data(MATRIX_RUN)
    on09 = surv[surv.survey == 2009]
    gis09 = on09.gis.to_numpy()
    idx09 = ((on09.gis - 1) * R.CELLS + on09.cell).to_numpy()
    def mod_line(k, idx=idx09, gis=gis09):
        return pd.Series(mod[k][idx]).groupby(gis).agg(stat)
    paired = {int(g.model_state.iloc[0][:4]) - first: (year, g) for year, g in surv.groupby("survey")}
    frames, last = [], None
    for k in range(mod.shape[0]):
        fresh = k in paired
        if fresh:
            year, g = paired[k]
            obs = g.groupby("gis").crest.agg(stat)
            m = mod_line(k, ((g.gis - 1) * R.CELLS + g.cell).to_numpy(), g.gis.to_numpy()).reindex(obs.index)
            d = (m - obs).dropna()
            last = (year, obs, f"r = {np.corrcoef(obs[d.index], m[d.index])[0, 1]:.2f}   "
                               f"RMSE = {np.sqrt((d ** 2).mean()):.2f} m   bias = {d.mean():+.2f} m")
        fig, ax = plt.subplots(figsize=fs.figsize("double", height=3.0), layout="constrained")
        line = mod_line(k)
        ax.plot(line.index, line, color=fs.INK, lw=1.6, zorder=4, label="Model")
        year, obs, note = last
        if fresh:
            ax.plot(obs.index, obs, color=C09, lw=1.0, marker="o", ms=3, zorder=5, label=f"Observed ({year} lidar)")
        else:
            ax.plot(obs.index, obs, color=C09, lw=0.8, alpha=0.35, zorder=3, label=f"Last survey ({year})")
        ax.set_title(f"1 January {first + k}", loc="left", fontsize=11, fontweight="bold", color=fs.INK)
        ax.set_title(note if fresh else " ", loc="right", fontsize=7, color=fs.INK_MUTED)
        ax.set(xlabel="Domain (GIS ID), south to north", ylabel=f"Domain-{stat} dune crest (m NAVD88)",
               xlim=(0.4, 90.6), ylim=(0, 10))
        ax.set_xticks(range(0, 91, 10))
        finish(ax)
        ax.legend(frameon=False, fontsize=7, loc="upper right", ncol=2)
        fig.set_dpi(GIF_DPI)
        fig.canvas.draw()
        frames.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba())).convert("RGB"))
        plt.close(fig)
    name = "stable_dune_fit_island_alongshore_1996_2009" + ("_mean" if stat == "mean" else "")
    path = OUT / "3-model-fit" / "island-wide" / f"{name}.gif"
    path.parent.mkdir(parents=True, exist_ok=True)
    ms = [int(GIF_FRAME_S * 1000)] * (len(frames) - 1) + [int(GIF_HOLD_S * 1000)]
    frames[0].save(path, save_all=True, append_images=frames[1:], duration=ms, loop=0, optimize=True)
    fs.record_caption(path, f"Animation of stable_dune_fit_island_alongshore through 1996-2009: domain-{stat} dune "
                      f"crest along Hatteras Island, one frame per model state (1 January 1997-2010, r 0.01-0.03, the "
                      f"current runner's 1996 full-management matrix run). Black: the model on the 10 m cells the 2009 "
                      f"lidar covers, the same cells in every frame. Amber: each lidar survey on the frame its scoring "
                      f"pairs it with (1 January the year after), kept faint until the next survey; the 1999 and 2000 "
                      f"surveys share the 2000 frame, where 2000 is drawn. r, RMSE and bias on that survey's own cells. Lidar spikes "
                      f"(crests more than 5 m above the same cell's before and after, up to 60 m, mostly 1996 and 1998) "
                      f"are left out upstream by HAT_survey_dune_crests.py. Surveys are irregular (none 2002-2003 or "
                      f"2006-2008). Producer: HAT_stable_dune_figures.py gif.")
    print(path, flush=True)


# Per rate and stretch, from the scored tables (HAT_stable_dune_growth_rate.py score, made with the runs):
# 2010 model minus 2009 lidar, and RMSE over the surveys after 1996; n = domains, for the pooled point
def fit_errors():
    t = pd.read_csv(R.OUT_DIR / "tables" / "stable_dune_summary.csv", header=[0, 1], index_col=0)
    rows = []
    for rate, row in t.iterrows():
        r = tuple(float(x) for x in rate.split("-"))
        for label, doms in STABLE.items():
            rows.append(dict(rate=r, group=label, n=len(doms), end=row[("2010 error vs 2009 lidar (m)", label)],
                             rmse=row[("RMSE over the surveys (m)", label)]))
    return pd.DataFrame(rows)


def fig_fit_errors():
    d = fit_errors()
    order = sorted(R.RATES, key=lambda r: (np.mean(r), r[1] - r[0]))
    y = {r: k for k, r in enumerate(order)}
    fig, axes = plt.subplots(1, 2, figsize=fs.figsize("double", height=3.4), layout="constrained", sharey=True,
                             gridspec_kw=dict(width_ratios=[1.5, 1]))
    for i, (ax, col, xl) in enumerate(zip(axes, ("end", "rmse"), ("Model minus lidar, median crest,\n2009 survey (m)",
                                                                  "RMSE of the median crest over\nthe surveys after 1996 (m)"))):
        ax.axhspan(y[ADOPTED] - 0.45, y[ADOPTED] + 0.45, color="0.9", zorder=0)
        for label in STABLE:
            g = d[d.group == label]
            ax.plot(g[col], g.rate.map(y), ls="none", marker="o", ms=4.5, color=STABLE_C[label],
                    label=f"Model minus lidar, {label}")
        pooled = d.groupby(d.rate.map(str)).apply(lambda g: np.average(g[col], weights=g.n))
        ax.plot([pooled[str(r)] for r in order], range(len(order)), ls="none", marker="D", ms=3.5, mfc="white",
                mec=fs.INK, mew=1.0, label="All seven domains")
        if col == "end":
            ax.axvline(0, color=fs.INK, lw=0.8, ls="--")
        ax.set(xlabel=xl)
        ax.grid(axis="x", color=fs.GRID_C, lw=0.5)
        fs.open_frame(ax)
        fs._letter_inside(ax, i)
    axes[0].set_yticks(range(len(order)), [rate_label(r) + ("  (adopted)" if r == ADOPTED else "") for r in order])
    axes[0].set_ylim(-0.6, len(order) - 0.4)
    axes[0].set_ylabel("Growth rate (rmin-rmax), by mean")
    fig.legend(*axes[0].get_legend_handles_labels(), frameon=False, fontsize=6.5, loc="outside lower center", ncol=3)
    save(fig, "stable_dune_fit_error_by_rate",
         "Every growth rate tested on the stable dunes, 1996 start (ALACE offset corrected), Dmax 10.5 m above the "
         "berm, groin on, ordered by mean rate. (a) Model minus lidar median crest at the 2009 survey, per stretch "
         "and for the seven domains together; 0 matches the lidar. (b) RMSE of the median crest over the "
         "surveys after 1996. Numbers from tables/stable_dune_summary.csv. Rates averaging 0.02 (0.02, 0.01-0.03, 0.005-0.035) are the ones near zero at both "
         "stretches with the lowest RMSE; slower rates leave the dunes too low, faster ones too high. Shaded: the "
         "adopted range.")


def fig_fit_spread():
    cells = stable_cells()
    obs, mods = [], {r: [] for r in SPREAD_RATES}
    runs = {r: model_crests(R.G.name(R.CEILING, r, groin=True)) for r in SPREAD_RATES}
    for label in STABLE:
        idx, s = cells[label]
        w = s.pivot_table(index=["gis", "cell"], values="crest", columns="survey")
        obs.append((w[2009] - w[1996]).to_numpy())
        st96 = s[s.survey == 1996].model_state.iloc[0]
        st09 = s[s.survey == 2009].model_state.iloc[0]
        for r, mod in runs.items():
            mods[r].append(model_at(mod, st09)[idx] - model_at(mod, st96)[idx])
    series = [("Observed: lidar, 1996 to 2009", np.concatenate(obs), "0.55")] + \
             [(lab, np.concatenate(mods[r]), fs.INK if r == ADOPTED else "#92c5de") for r, lab in SPREAD_RATES.items()]
    fig, ax = plt.subplots(figsize=fs.figsize("double", height=3.0), layout="constrained")
    for k, (lab, v, col) in enumerate(series):
        y = len(series) - 1 - k
        parts = ax.violinplot(v, positions=[y], vert=False, widths=0.8, showextrema=False)
        for b in parts["bodies"]:
            b.set_facecolor(col)
            b.set_alpha(0.35)
            b.set_edgecolor("none")
        q1, med, q3 = np.percentile(v, [25, 50, 75])
        ax.plot([q1, q3], [y, y], color=col, lw=3, solid_capstyle="butt")
        ax.plot(med, y, marker="o", ms=4.5, mfc="white", mec=col, mew=1.3)
    ax.axvline(np.median(series[0][1]), color="0.55", lw=0.6, ls="--")
    ticks = [f"{lab}\nmedian {np.median(v):+.2f} m, middle half {np.subtract(*np.percentile(v, [75, 25])):.2f} m wide"
             for lab, v, _ in series]
    ax.set_yticks(range(len(series)), ticks[::-1], fontsize=6.5)
    ax.set(xlabel="Crest change per 10 m cell, 1996-2009 (m)", xlim=(-1.0, 2.5))
    ax.grid(axis="x", color=fs.GRID_C, lw=0.5)
    fs.open_frame(ax)
    save(fig, "stable_dune_fit_spread",
         f"Why a range of growth rates rather than one: crest change per 10 m cell at the seven stable domains, "
         f"1996 to 2009 ({len(series[0][1])} cells with a crest in every survey). Shape: distribution; bar: middle "
         f"half; circle: median; dashed: the observed median. All three model rates match the observed median, but a "
         f"fixed 0.02 grows every cell alike. 0.01-0.03 doubles the cell-to-cell spread at the same median; "
         f"0.005-0.035 widens it further but lowers the median. Barrier3D draws each cell's r at random within the "
         f"range, so the spread is not placed where the real fast and slow dunes are.")


def main():
    fs.apply_style()
    OUT.mkdir(parents=True, exist_ok=True)
    if len(sys.argv) > 1 and sys.argv[1] == "fit":            # only the why-this-rate figures
        fig_fit_trajectory()
        fig_fit_errors()
        fig_fit_spread()
        fig_fit_island()
        fig_fit_island_alongshore()
        return
    if len(sys.argv) > 1 and sys.argv[1] == "gif":            # only the alongshore animation, 1996-2009
        gif_fit_island_alongshore(stat="mean")
        return
    if len(sys.argv) > 1 and sys.argv[1] == "alongshore":     # only the island-wide alongshore figure
        fig_fit_island_alongshore()
        return
    if len(sys.argv) > 1 and sys.argv[1] == "dems":           # only the same-dune and per-domain figures (no model runs)
        for g in sorted(g for doms in STABLE.values() for g in doms):
            fig_domain(g)
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
    for g in sorted(g for doms in STABLE.values() for g in doms):
        fig_domain(g)
    fig_alongshore()
    fig_crest_shift()
    fig_fit_trajectory()
    fig_fit_errors()
    fig_fit_spread()
    fig_fit_island()
    fig_fit_island_alongshore()


if __name__ == "__main__":
    main()
