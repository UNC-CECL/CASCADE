#!/usr/bin/env python3
"""
House-style figures for the shoreline rates around the Buxton groins, from the tables HAT_groin_shoreline_analysis.py writes.

    python shoreline_rates_by_era_figures.py   ->  figures/shoreline_rates_by_era_<n>_<what>.png

Four figures, one question each: how the rate differed in each era of the groin's life,
what the groin changed and how far, how the rate evolved decade by decade, and how each
side's mean rate near the groin changed. Each panel title is the question; captions (what
it tests, how to read it, what it shows) go to figures/supporting/CAPTIONS.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-08
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import BoundaryNorm  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from statsmodels.nonparametric.smoothers_lowess import lowess  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = next(p for p in HERE.parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))
from site_layer.hat_figure_style import (  # noqa: E402
    C, C_1984, C_1997, DOMAIN_AXIS_LABEL, INK, INK_MUTED, _letter_inside, _title, apply_style,
    figsize, open_frame, record_caption, save, town_bands)

# --- CONFIG ------------------------------------------------------------------
TABLES = HERE / "output"
FIGS = HERE / "figures"
PRE, WORK, AFTER = "Pre-install", "Functional groin", "Deteriorated"
ERA = {PRE: ("before the groin, 1849-1969", C["BASE"], "--"),
       WORK: ("groin working, 1970-1995", C_1984, "-"),
       AFTER: ("after the 1995 repair, 1996-2024", C_1997, "-")}
SMOOTH_DOMAINS = 7.0                    # the group's LOWESS width for observations, whole reach
NEAR_SMOOTH_DOMAINS = 2.0               # near the groins: the south side is only ~5.5 domains long
GAP_DOMAINS = 1.5                       # no smoothed line across a stretch without transects
ZONES = {"downdrift": (1, 4), "updrift": (7, 20)}   # as the analysis script's zone panels
NEAR = (1, 20)                          # the near-groin view
RATE_LIM = 5.0                          # colour limit of the decade map, m/yr
EXTENT_DECADES = (1970, 1980, 1990, 2000, 2010, 2020)
# -----------------------------------------------------------------------------


def tab(name):
    return pd.read_csv(TABLES / f"groin_analysis_{name}.csv")


# Distance from the groin (m) <-> a continuous GIS domain coordinate, from per-domain means
def domain_axis(t):
    agg = t.groupby("domain")["dist_from_groin_m"].mean().sort_index()
    d, m = agg.index.values.astype(float), agg.values
    order = np.argsort(m)

    def to_dom(dist_m):
        return np.interp(dist_m, m[order], d[order])

    def to_km(dom):
        return np.interp(dom, d, m) / 1000.0
    return to_dom, to_km


def _lowess(x, y, width):
    if len(x) < 10:
        return np.array([]), np.array([])
    span = max(x.max() - x.min(), 1e-6)
    frac = min(1.0, min(width, span / 3) / span)       # a short side gets a third of its length
    s = lowess(y, x, frac=frac, return_sorted=True)
    xs, ys = s[:, 0], s[:, 1]
    gap = np.r_[False, np.diff(xs) > GAP_DOMAINS]
    return np.insert(xs, np.where(gap)[0], np.nan), np.insert(ys, np.where(gap)[0], np.nan)


# LOWESS on each side of the groins separately, so the two sides are never blended
def smooth(x, y, gx, width=SMOOTH_DOMAINS):
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    parts = [_lowess(x[m], y[m], width) for m in (x < gx, x >= gx)]
    return np.r_[parts[0][0], np.nan, parts[1][0]], np.r_[parts[0][1], np.nan, parts[1][1]]


def groin_line(ax, x):
    ax.axvline(x, color=C["GROIN"], lw=1.0, zorder=4)


def km_axis(ax, to_dom, to_km):
    sec = ax.secondary_xaxis("top", functions=(to_km, lambda k: to_dom(np.asarray(k) * 1000.0)))
    lo, hi = sorted(to_km(np.array(ax.get_xlim())))
    step = next(st for st in (1, 2, 5, 10, 20) if (hi - lo) / st <= 8)
    sec.set_xticks(np.arange(np.ceil(lo / step) * step, hi + 1e-9, step))
    sec.set_xlabel("Distance from the groins (km, + = north, updrift)", fontsize=7.5)
    sec.tick_params(labelsize=7)


def alongshore(ax, lim, to_dom, to_km, km=True):
    ax.set_xlim(*lim)
    ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    town_bands(ax, strip=0.05)
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.grid(axis="y")
    open_frame(ax)
    if km:
        km_axis(ax, to_dom, to_km)


# 1. The rate in each era along the coast
def fig_eras(e, to_dom, to_km, gx):
    fig, axes = plt.subplots(2, 1, figsize=figsize("double", height=6.4), constrained_layout=True)
    full = (np.floor(e["x"].min()), np.ceil(e["x"].max()))
    for i, (ax, lim, title) in enumerate(zip(axes, (full, (NEAR[0] - 0.5, NEAR[1] + 0.5)), (
            "How did the rate differ before, during and after the groin's working life?",
            f"Close to the groins, GIS {NEAR[0]}-{NEAR[1]}"))):
        for era, (lab, col, ls) in ERA.items():
            s = e[(e["era"] == era) & e["x"].between(*lim)]
            ax.plot(s["x"], s["slope_m_yr"], ls="none", marker=".", ms=2.0, color=col, alpha=0.25, zorder=2)
            xs, ys = smooth(s["x"].to_numpy(), s["slope_m_yr"].to_numpy(), gx,
                            SMOOTH_DOMAINS if i == 0 else NEAR_SMOOTH_DOMAINS)
            ax.plot(xs, ys, color=col, ls=ls, lw=1.6, zorder=3)
        groin_line(ax, gx)
        alongshore(ax, lim, to_dom, to_km)
        ax.set_ylabel("Shoreline change rate (m/yr)\n+ = accretion")
        v = e[e["x"].between(*lim)]["slope_m_yr"]
        lo_v, hi_v = np.nanpercentile(v, [1, 99])
        pad = 0.15 * (hi_v - lo_v)
        ax.set_ylim(lo_v - pad, hi_v + pad)
        _title(ax, i, title)
    axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)
    handles = [Line2D([], [], color=col, ls=ls, lw=1.6, label=lab) for lab, col, ls in ERA.values()]
    handles += [Line2D([], [], ls="none", marker=".", ms=5, color=C["BASE"], alpha=0.5,
                       label="single transects (faint)"),
                Line2D([], [], color=C["GROIN"], lw=1.0, label="Buxton groins")]
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False)
    png = FIGS / "shoreline_rates_by_era_1_rates_by_era.png"
    save(fig, png, close=True)
    n = e.groupby("era")["transect_id"].nunique()
    record_caption(png, (
        "**Tests:** whether the shoreline change rate near Buxton differed between the decades before "
        "the groins were built, their working life, and the years after the 1995 last repair, and over "
        "what stretch of coast. **How to read:** each dot is the linear-regression rate (LRR) of one "
        "CoastSat transect over one era, positive for accretion; lines are LOWESS fits over "
        f"{SMOOTH_DOMAINS:g} domains in (a) and {NEAR_SMOOTH_DOMAINS:g} in (b), fitted separately on each side "
        "of the groins and broken where there are no transects; south of the groins there are only about "
        "5.5 domains of coast, too short for the 7-domain window, so on any side the window is capped at a third of its length. Eras: before the groin "
        f"(1849-1969, {n.get(PRE, 0)} transects with historical shorelines), the groin working "
        f"(1970-1995, {n.get(WORK, 0)}), after the last repair (1996-2024, {n.get(AFTER, 0)}). The "
        "pre-install rates rest on three historical surveys per transect (1849/1852/1860 or "
        "1946/1949/1967), so their fits are near-perfect by construction and say little about their "
        "uncertainty. Red line: the groin field; grey strips: villages. (a) The whole analysed reach. "
        f"(b) GIS {NEAR[0]}-{NEAR[1]}. The top axis gives distance from the groins."))


# 2. What the groin changed, era to era
def fig_differences(e, to_dom, to_km, gx):
    w = e.pivot_table(index=["transect_id"], columns="era", values="slope_m_yr")
    x = e.groupby("transect_id")["x"].first().reindex(w.index)
    pairs = ((WORK, PRE, "What did the groin's working life change?",
              "Working life minus before the groin (m/yr)", C_1984),
             (AFTER, WORK, "What changed once the groin stopped holding?",
              "After the repair minus working life (m/yr)", C_1997))
    fig, axes = plt.subplots(2, 1, figsize=figsize("double", height=6.0), sharex=True, constrained_layout=True)
    lim = (NEAR[0] - 0.5, 30.5)
    stats = []
    for i, (ax, (a, b, title, ylab, col)) in enumerate(zip(axes, pairs)):
        d = (w[a] - w[b]).dropna()
        xx = x.reindex(d.index)
        keep = xx.between(*lim)
        ax.plot(xx[keep], d[keep], ls="none", marker=".", ms=2.2, color=col, alpha=0.3)
        xs, ys = smooth(xx[keep].to_numpy(), d[keep].to_numpy(), gx, NEAR_SMOOTH_DOMAINS)
        ax.plot(xs, ys, color=col, lw=1.6)
        groin_line(ax, gx)
        alongshore(ax, lim, to_dom, to_km, km=(i == 0))
        ax.set_ylabel(ylab.replace(" (m/yr)", "\n(m/yr)"))
        _title(ax, i, title)
        for side, (lo, hi) in ZONES.items():
            z = d[xx.between(lo - 0.5, hi + 0.5)]
            stats.append((a, b, side, z.mean(), len(z)))
    axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)
    handles = [Line2D([], [], ls="none", marker=".", ms=5, color=INK_MUTED, label="single transects (faint)"),
               Line2D([], [], color=INK_MUTED, lw=1.6, label=f"LOWESS over {NEAR_SMOOTH_DOMAINS:g} domains, each side"),
               Line2D([], [], color=C["GROIN"], lw=1.0, label="Buxton groins")]
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False)
    png = FIGS / "shoreline_rates_by_era_2_era_differences.png"
    save(fig, png, close=True)
    s = {(a, side): (m, k) for a, b, side, m, k in stats}
    record_caption(png, (
        "**Tests:** how much the groin changed the shoreline change rate on each side, and whether that "
        "change reversed once it stopped holding. **How to read:** per transect, the rate in one era "
        "minus the rate in the era before it, positive where the shoreline accreted faster (or eroded "
        "more slowly) in the later era; line: LOWESS over "
        f"{NEAR_SMOOTH_DOMAINS:g} domains, fitted separately on each side of the groins. A working groin "
        "should raise the rate on the north (updrift) side "
        "and lower it on the south (downdrift) side. (a) Working life (1970-1995) minus before the "
        "groin (1849-1969). (b) After the repair (1996-2024) minus working life. Red line: the groin "
        f"field. **Shows:** mean difference in (a) {s[(WORK, 'downdrift')][0]:+.1f} m/yr south of the "
        f"groins (GIS {ZONES['downdrift'][0]}-{ZONES['downdrift'][1]}, {s[(WORK, 'downdrift')][1]} "
        f"transects) and {s[(WORK, 'updrift')][0]:+.1f} m/yr north (GIS {ZONES['updrift'][0]}-"
        f"{ZONES['updrift'][1]}, {s[(WORK, 'updrift')][1]}); in (b) {s[(AFTER, 'downdrift')][0]:+.1f} "
        f"and {s[(AFTER, 'updrift')][0]:+.1f} m/yr. Because the pre-install rates come from three "
        "historical surveys, (a) carries more uncertainty than (b)."))


# 3. Decade by decade, with how far the groin's effect reached
def fig_decades(dec, ext, to_dom, to_km, gx):
    dec = dec[dec["decade_start"] >= 1970].copy()
    decades = sorted(dec["decade_start"].unique())
    lim = (NEAR[0] - 0.5, 40.5)
    fig, ax = plt.subplots(figsize=figsize("double", height=3.8), constrained_layout=True)
    levels = np.arange(-RATE_LIM, RATE_LIM + 0.5, 1.0)
    cmap = plt.get_cmap("RdBu", len(levels) - 1)
    norm = BoundaryNorm(levels, cmap.N)
    for k, dstart in enumerate(decades):
        s = dec[(dec["decade_start"] == dstart) & dec["x"].between(*lim)].sort_values("x")
        if s.empty:
            continue
        xe = np.r_[s["x"].iloc[0] - 0.05, (s["x"].values[1:] + s["x"].values[:-1]) / 2, s["x"].iloc[-1] + 0.05]
        xe = np.clip(xe, np.r_[xe[0], s["x"].values - 0.15], np.r_[s["x"].values + 0.15, xe[-1]])
        ax.pcolormesh(xe, [k - 0.45, k + 0.45], s["slope_m_yr"].values[None, :], cmap=cmap, norm=norm,
                      shading="flat")
    for k, dstart in enumerate(decades):
        r = ext[ext["decade_start"] == dstart]
        if r.empty:
            continue
        r = r.iloc[0]
        for col_m, mk in (("updrift_extent_m", ">"), ("downdrift_extent_m", "<")):
            if np.isfinite(r[col_m]) and r[col_m] > 0:
                sign = 1 if col_m.startswith("up") else -1
                ax.plot(float(to_dom(sign * r[col_m])), k, marker=mk, ms=5, color=INK, zorder=6)
    groin_line(ax, gx)
    ax.set_xlim(*lim)
    ax.set_ylim(-1.0, len(decades) - 0.4)
    ax.set_yticks(range(len(decades)), [f"{d}s" for d in decades])
    ax.invert_yaxis()
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    ax.set_ylabel("Decade")
    town_bands(ax, strip=0.06)
    ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    km_axis(ax, to_dom, to_km)
    open_frame(ax)
    ax.set_title("How did the rate evolve decade by decade, and how far did the groin's effect reach?")
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, pad=0.01, extend="both")
    cb.set_label("Shoreline change rate (m/yr)\n+ = accretion")
    handles = [Line2D([], [], ls="none", marker=">", ms=5, color=INK, label="end of the groin's effect, north"),
               Line2D([], [], ls="none", marker="<", ms=5, color=INK, label="end of the groin's effect, south"),
               Line2D([], [], color=C["GROIN"], lw=1.0, label="Buxton groins")]
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False)
    png = FIGS / "shoreline_rates_by_era_3_decade_map.png"
    save(fig, png, close=True)
    e70 = ext[ext["decade_start"] == 1970].iloc[0]
    record_caption(png, (
        "**Tests:** when, decade by decade, the shoreline near the groins accreted or eroded, and how far "
        "along the coast the groin's effect reached in each decade. **How to read:** one row per decade; "
        "each column is one CoastSat transect, coloured by its linear-regression rate over that decade "
        f"(red erosion, blue accretion, white near zero; capped at ±{RATE_LIM:g} m/yr); blank where a "
        "transect has too few images in that decade. Triangles: where the decade's rate stops differing "
        "from the pre-install baseline by more than 1 m/yr, searched outward from the groins in 500 m "
        "bins (the analysis script's signal extent). Red line: the groin field. **Shows:** in the 1970s "
        f"the effect reached {e70.updrift_extent_m / 1000:.1f} km north and "
        f"{e70.downdrift_extent_m / 1000:.1f} km south of the groins; the 1970s row rests on few images "
        "(Landsat coverage begins in the mid-1980s for most transects)."))


# 4. Zone means near the groin, era by era
def fig_zones(e, inc5, inc10):
    rows = [(PRE, e, PRE), ("first 5 years", inc5, "1970-1974"), ("first 10 years", inc10, "1970-1979"),
            (WORK, e, WORK), (AFTER, e, AFTER)]
    labels = {PRE: "before the groin\n1849-1969", "first 5 years": "first 5 years\n1970-1974",
              "first 10 years": "first 10 years\n1970-1979", WORK: "working life\n1970-1995",
              AFTER: "after the repair\n1996-2024"}
    out = []
    for key, t, era in rows:
        s = t[t["era"] == era]
        for side, (lo, hi) in ZONES.items():
            z = s[s["x"].between(lo - 0.5, hi + 0.5)]["slope_m_yr"].dropna()
            out.append(dict(period=key, side=side, mean=z.mean(), se=z.std(ddof=1) / np.sqrt(len(z)) if len(z) > 1
                            else np.nan, n=len(z)))
    o = pd.DataFrame(out)
    fig, ax = plt.subplots(figsize=figsize("double", height=3.2), constrained_layout=True)
    keys = [r[0] for r in rows]
    side_col = {"downdrift": C["ADDED"], "updrift": C_1997}
    for k, side in enumerate(("downdrift", "updrift")):
        s = o[o["side"] == side].set_index("period").reindex(keys)
        x = np.arange(len(keys)) + (k - 0.5) * 0.24
        ax.errorbar(x, s["mean"], yerr=s["se"], fmt="o", ms=5, capsize=2.5, lw=1.0, color=side_col[side],
                    label=f"{'south (downdrift)' if side == 'downdrift' else 'north (updrift)'}, "
                          f"GIS {ZONES[side][0]}-{ZONES[side][1]}")
        for xi, (_, r) in zip(x, s.iterrows()):
            ax.text(xi, r["mean"] + (r["se"] if np.isfinite(r["se"]) else 0) + 0.3, f"n={int(r['n'])}",
                    ha="center", fontsize=6, color=side_col[side])
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.axvspan(0.5, 3.5, color=C_1984, alpha=0.05, lw=0)
    ax.set_xticks(range(len(keys)), [labels[k] for k in keys])
    ax.set_ylabel("Mean shoreline change rate\n(m/yr, + = accretion, ± 1 SE)")
    ax.grid(axis="y")
    open_frame(ax)
    ax.set_title("Near the groins, how did each side's mean rate change?")
    fig.legend(loc="outside lower center", ncol=2, frameon=False)
    png = FIGS / "shoreline_rates_by_era_4_zone_means.png"
    save(fig, png, close=True)
    o.to_csv(FIGS / "supporting" / "zone_means.csv", index=False, float_format="%.3f")
    g = o.set_index(["period", "side"])
    record_caption(png, (
        "**Tests:** whether the groin made the north (updrift) side accrete and the south (downdrift) "
        "side erode relative to before it was built, how quickly, and whether that lasted. **How to "
        f"read:** the mean linear-regression rate of the transects in two zones, GIS {ZONES['downdrift'][0]}-"
        f"{ZONES['downdrift'][1]} south of the groins and GIS {ZONES['updrift'][0]}-{ZONES['updrift'][1]} "
        "north, ± 1 standard error across transects; n is the number of transects with a rate in that "
        "period. The tinted band is the groin's working life; the first 5 and 10 years are its opening "
        "windows, fitted on far fewer images. **Shows:** north of the groins the mean rate goes from "
        f"{g.loc[(PRE, 'updrift'), 'mean']:+.1f} m/yr before the groin to "
        f"{g.loc[('first 5 years', 'updrift'), 'mean']:+.1f} in its first 5 years and "
        f"{g.loc[(WORK, 'updrift'), 'mean']:+.1f} over its working life, then "
        f"{g.loc[(AFTER, 'updrift'), 'mean']:+.1f} after the repair; south of the groins "
        f"{g.loc[(PRE, 'downdrift'), 'mean']:+.1f}, {g.loc[('first 5 years', 'downdrift'), 'mean']:+.1f}, "
        f"{g.loc[(WORK, 'downdrift'), 'mean']:+.1f} and {g.loc[(AFTER, 'downdrift'), 'mean']:+.1f} m/yr. "
        "Table: supporting/zone_means.csv."))


def main():
    apply_style()
    (FIGS / "supporting").mkdir(parents=True, exist_ok=True)
    e = tab("era_lrrs")
    to_dom, to_km = domain_axis(e)
    gx = float(to_dom(0.0))
    pos = e.groupby("transect_id")["dist_from_groin_m"].first()
    e["x"] = to_dom(e["dist_from_groin_m"])
    inc5, inc10 = tab("decade_increment_lrrs_5yr"), tab("decade_increment_lrrs_10yr")
    for t in (inc5, inc10):
        t["x"] = to_dom(t["dist_from_groin_m"])
    dec = tab("decadal_lrr")
    dec["x"] = to_dom(dec["transect_id"].map(pos))
    ext = tab("signal_extent")
    fig_eras(e, to_dom, to_km, gx)
    fig_differences(e, to_dom, to_km, gx)
    fig_decades(dec, ext, to_dom, to_km, gx)
    fig_zones(e, inc5, inc10)
    print("drew 4 figures")


if __name__ == "__main__":
    main()
