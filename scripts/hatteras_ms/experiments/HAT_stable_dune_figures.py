"""
Why the stable dunes at GIS 59-62 and 20-22 can calibrate the dune growth rate r, one point per figure.

    python scripts/hatteras_ms/experiments/HAT_stable_dune_figures.py

Figures for the stable-dune study (HAT_stable_dune_growth_rate.py, run it first),
written to its figures/ folder (data/.../1996-start/dune-growth-rate/ since 2026-10-21) in reading order:
  stable_dune_fig1_where       shoreline change, dune change and overwash per domain, 1996-2009:
                               the two stretches meet all three conditions
  stable_dune_fig2_profiles    1996 and 2009 DEM profiles at the stable stretches: dune and beach in place
  stable_dune_fig3_excluded    the same at two excluded domains: a growing beach, a retreating shoreline
  stable_dune_fig4_surveys     crests through all nine 1996-2009 lidar surveys: no storm losses
  stable_dune_fig5_model       model crest for several r against the surveys
  stable_dune_fig6_rate        (drawn by HAT_stable_dune_growth_rate.py) model minus lidar against r
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
    path = OUT / f"{name}.png"
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


def main():
    fs.apply_style()
    OUT.mkdir(parents=True, exist_ok=True)
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
