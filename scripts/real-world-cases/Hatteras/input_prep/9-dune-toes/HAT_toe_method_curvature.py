"""
Maximum curvature (Stockdon, Sallenger, Holman & Howd 2007, Marine Geology 238, 1-20, "A simple model for the spatially-
variable coastal response to hurricanes"; as implemented by pybeach, Beuzen 2019): the dune toe is where the slope
increases fastest landward, the largest second derivative d2z/dx2 of the profile between the shoreline and the dune
crest. Here on our profiles, shoreline and crest (HAT_toe_method_line.run), the profile smoothed over 5 m, searched
in the lower half of the rise (unconstrained, a kink on the dune face wins), no lower than the berm as ours.

    python scripts/real-world-cases/Hatteras/input_prep/9-dune-toes/HAT_toe_method_curvature.py                # the method figure (stockdon_*.png)
    python scripts/real-world-cases/Hatteras/input_prep/9-dune-toes/HAT_toe_method_curvature.py --check        # against the line method, three smoothings
    python scripts/real-world-cases/Hatteras/input_prep/9-dune-toes/HAT_toe_method_curvature.py --why          # why it fails on Hatteras, one profile

Writes, in data/hatteras_init/9-dune-toes/1-transects-1996/3-methods/3-curvature/ (D:): stockdon_a-d.png and
stockdon_method_1996.png (the method figure), curvature_check_1996.csv (domains) and
curvature_check_transects_1996.csv with curvature_check.png, curvature_why.png.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-05
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parent))
from HAT_toe_common import (BERM_M, CHORD, CREST, EXAMPLE_GIS, FILL, MHW, MIN_KEPT, M_CURVATURE, PROFILE, STEP_M,
    TOE, VIEW_M, WATER, YEAR, _legend)  # noqa: E402
from HAT_toe_method_line import (run)  # noqa: E402

CURV_SMOOTH_M = (3, 5, 9)        # running-mean widths tried before the second derivative; 5 m is the main one
CURV_MAIN_M = 5


def curvature_toe(s, z, p, smooth_m=CURV_MAIN_M, lower=True):
    """The maximum-curvature toe (the published alternative; Stockdon et al. 2007 as cited for it): on the same
    profile, shoreline and crest as pick(), the point between them where the slope increases fastest landward
    (largest d2z/dx2 of the profile smoothed by a smooth_m running mean), no lower than the berm as pick().
    lower: only where the profile is in the lower half between the shoreline and the crest (a toe is not up the dune
    face; unconstrained, a kink or swale on the face is the sharpest bend and wins)."""
    if p.get("toe") is None or p.get("crest") is None:
        return None
    i0, i1 = p["shore"], p["crest"]
    k = max(1, int(round(smooth_m / STEP_M)))
    zz = pd.Series(z[i0:i1 + 1]).rolling(k, center=True, min_periods=1).mean().to_numpy()
    if len(zz) < 5 or not np.isfinite(zz).all():
        return None
    d2 = np.gradient(np.gradient(zz, STEP_M), STEP_M)
    d2[:2], d2[-2:] = -np.inf, -np.inf                           # not the shoreline or crest themselves
    if lower:
        d2[zz > zz[0] + 0.5 * (zz[-1] - zz[0])] = -np.inf
    toe = i0 + int(np.argmax(d2))
    if z[toe] < BERM_M:
        toe = next((j for j in range(toe, i1 + 1) if z[j] >= BERM_M), i1)
    return toe


def curvature_check(t, dom, profiles, *_):
    """The maximum-curvature toe against the line method, transect by transect and as domain medians, at three
    smoothing widths. Writes curvature_check_<year>.csv (domains) and figures/curvature_check.png."""
    rows = []
    for r in t.itertuples():
        s, z, p = profiles[r.transect][:3]
        row = dict(transect=r.transect, gis=r.gis, line_toe=r.toe_m_navd88 if r.kept else np.nan)
        for w in CURV_SMOOTH_M:
            i = curvature_toe(s, z, p, w) if r.kept else None
            row[f"curv{w}_toe"] = float(z[i]) if i is not None else np.nan
            row[f"curv{w}_x"] = float(s[i] - s[p["shore"]]) if i is not None else np.nan
        i = curvature_toe(s, z, p, CURV_MAIN_M, lower=False) if r.kept else None
        row["curvall_toe"] = float(z[i]) if i is not None else np.nan
        row["line_x"] = float(r.toe_from_shoreline_m) if r.kept else np.nan
        rows.append(row)
    c = pd.DataFrame(rows)
    med = lambda col: c.groupby("gis")[col].median().where(c.groupby("gis")[col].count() >= MIN_KEPT)
    d = pd.DataFrame({"line": med("line_toe"), "curvall": med("curvall_toe")}
                     | {f"curv{w}": med(f"curv{w}_toe") for w in CURV_SMOOTH_M})
    d.index.name = "gis"
    (M_CURVATURE).mkdir(parents=True, exist_ok=True)
    d.round(3).to_csv(M_CURVATURE / f"curvature_check_{YEAR}.csv")
    c.round(3).to_csv(M_CURVATURE / f"curvature_check_transects_{YEAR}.csv", index=False)
    for w in ("all",) + CURV_SMOOTH_M:
        dt = (c[f"curv{w}_toe"] - c.line_toe).dropna()
        dd = (d[f"curv{w}"] - d.line).dropna()
        print(f"curvature, {'unconstrained, 5' if w == 'all' else str(w) + ' m, lower half,'} m smoothing: transects within 0.25 m of the line toe {100 * (dt.abs() <= 0.25).mean():.0f}%,"
              f" within 0.5 m {100 * (dt.abs() <= 0.5).mean():.0f}% (n {len(dt)}, mean {dt.mean():+.2f} m); "
              f"domain medians: mean {dd.mean():+.2f} m, RMSE {np.sqrt((dd ** 2).mean()):.2f} m, "
              f"> 0.5 m apart in {int((dd.abs() > 0.5).sum())} ({', '.join(str(g) for g in dd.index[dd.abs() > 0.5])})")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer import hat_figure_style as fs
    fs.apply_style()
    fig, axes = plt.subplots(2, 1, figsize=(fs.figsize("double")[0], 5.6), layout="constrained")
    # (a) the example transect: both picks, and the curvature that sets the second
    ax = axes[0]
    ids = list(t[t.gis == EXAMPLE_GIS].transect)
    ex = ids[len(ids) // 2]
    s, z, p = profiles[ex][:3]
    x = s - s[p["shore"]]
    v = (x >= VIEW_M[0]) & (x <= VIEW_M[1])
    ax.fill_between(x[v], z[v], -1, color=FILL, lw=0)
    ax.plot(x[v], z[v], color=PROFILE, lw=1.1, label=f"Profile, transect {ex}")
    ax.axhline(MHW, color=WATER, lw=0.7)
    ax.axhline(BERM_M, color=CREST, lw=0.6, ls=(0, (2, 2)))
    ax.text(VIEW_M[0] + 2, BERM_M + 0.12, f"Berm ({BERM_M:g} m)", fontsize=7, color=CREST)
    ax.plot(x[p["crest"]], z[p["crest"]], "^", ms=7, color=CREST, label="Foredune crest")
    ax.plot(x[p["toe"]], z[p["toe"]], "o", ms=8, color=TOE, label=f"Line method ({z[p['toe']]:.2f} m)")
    iu = curvature_toe(s, z, p, lower=False)
    if iu is not None:
        ax.plot(x[iu], z[iu], "D", ms=6, mfc="white", mec="#E07B39", mew=1.0, zorder=6,
                label=f"Maximum curvature, whole face ({z[iu]:.2f} m)")
    ic = curvature_toe(s, z, p)
    if ic is not None:
        ax.plot(x[ic], z[ic], "D", ms=6, color="#E07B39", mec="white", mew=0.6, zorder=6,
                label=f"Maximum curvature, lower half ({z[ic]:.2f} m)")
    i0, i1 = p["shore"], p["crest"]
    zz = pd.Series(z[i0:i1 + 1]).rolling(CURV_MAIN_M, center=True, min_periods=1).mean().to_numpy()
    d2 = np.gradient(np.gradient(zz))
    ax2 = ax.twinx()
    ax2.plot(x[i0:i1 + 1], d2, color="#E07B39", lw=0.8, alpha=0.6)
    ax2.set_ylabel("Curvature d$^2$z/dx$^2$ (m$^{-1}$)", color="#E07B39", fontsize=8)
    ax2.tick_params(axis="y", colors="#E07B39", labelsize=7)
    ax2.spines["top"].set_visible(False)
    ax.set_xlim(*VIEW_M)
    ax.set_ylim(-0.5, np.ceil(np.nanmax(z[v]) + 0.8))
    ax.set_xlabel("Distance landward of the shoreline (m)")
    ax.set_ylabel("Elevation (m NAVD88)")
    _legend(ax)
    ax.text(-0.02, 1.02, "a", transform=ax.transAxes, fontsize=10, fontweight="bold", ha="right", va="bottom")
    # (b) domain medians alongshore
    ax = axes[1]
    ax.plot(d.index, d.line, color=TOE, lw=1.6, label="Line method (used)")
    ax.plot(d.index, d[f"curv{CURV_MAIN_M}"], color="#E07B39", lw=1.3,
            label=f"Maximum curvature, lower half ({CURV_MAIN_M} m smoothing)")
    for w, ls in zip([w for w in CURV_SMOOTH_M if w != CURV_MAIN_M], ((0, (1.5, 1.5)), (0, (4, 2)))):
        ax.plot(d.index, d[f"curv{w}"], color="#F4B183", lw=0.9, ls=ls, label=f"Lower half, {w} m smoothing")
    ax.set_xlim(0.5, 90.5)
    ax.set_xlabel(getattr(fs, "DOMAIN_AXIS_LABEL", "GIS domain (south → north)"))
    ax.set_ylabel("Dune toe (m NAVD88)")
    ax.set_ylim(np.floor(d.min().min() - 0.3), np.ceil(d.max().max() + 1.0))
    fs.open_frame(ax)
    _legend(ax, ncol=2)
    ax.text(-0.02, 1.02, "b", transform=ax.transAxes, fontsize=10, fontweight="bold", ha="right", va="bottom")
    png = M_CURVATURE / "curvature_check.png"
    fs.save(fig, png, close=True)
    dd = (d[f"curv{CURV_MAIN_M}"] - d.line).dropna()
    fs.record_caption(png, (
        f"The dune toe by maximum curvature against the line method, {YEAR} lidar. (a) Transect {ex} (GIS "
        f"{EXAMPLE_GIS}): the line method's toe (the point furthest below the shoreline-to-crest line) and the "
        f"maximum-curvature toe (the largest d2z/dx2, profile smoothed over {CURV_MAIN_M} m; orange line, right "
        f"axis), searched over the whole shoreline-to-crest face (open diamond) and over its lower half (filled), "
        f"all no lower than the berm. (b) Domain medians alongshore, the "
        f"curvature toe at {', '.join(str(w) for w in CURV_SMOOTH_M)} m smoothing. Curvature minus line: mean "
        f"{dd.mean():+.2f} m, RMSE {np.sqrt((dd ** 2).mean()):.2f} m."))
    print(f"wrote {png}")


def stockdon_figures(t, dom, profiles, *_, gis=EXAMPLE_GIS):
    """The method figure redone with the maximum-curvature toe (Stockdon et al. 2007 as cited for it; the whole
    shoreline-to-crest face, CURV_MAIN_M smoothing, berm floor), beside the line method: (a) one transect with its
    curvature, (b) the domain's five transects, (c) the domain toe, (d) domain toes alongshore.
    Writes figures/stockdon_method_<year>.png and one file per panel (stockdon_a ... d)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer import hat_figure_style as fs
    fs.apply_style()
    ORANGE = "#E07B39"
    pick_c = {r.transect: (curvature_toe(*profiles[r.transect][:3], lower=False) if r.kept else None)
              for r in t.itertuples()}
    ct = pd.Series({i: float(profiles[i][1][k]) for i, k in pick_c.items() if k is not None})
    t = t.assign(curv_toe_m_navd88=t.transect.map(ct))
    cnt = t.groupby("gis").curv_toe_m_navd88.count()
    cmed = t.groupby("gis").curv_toe_m_navd88.median().where(cnt >= MIN_KEPT)
    ids = list(t[t.gis == gis].transect)
    ex = ids[len(ids) // 2]
    shades = plt.get_cmap("Oranges")(np.linspace(0.4, 0.9, len(ids)))

    def pa(ax):                                                  # one transect and its curvature
        s, z, p = profiles[ex][:3]
        x = s - s[p["shore"]]
        v = (x >= VIEW_M[0]) & (x <= VIEW_M[1])
        i0, i1 = p["shore"], p["crest"]
        zz = pd.Series(z[i0:i1 + 1]).rolling(CURV_MAIN_M, center=True, min_periods=1).mean().to_numpy()
        d2 = np.gradient(np.gradient(zz, STEP_M), STEP_M)
        ax2 = ax.twinx()
        ax2.axhline(0, color=ORANGE, lw=0.4, alpha=0.5)
        ax2.plot(x[i0:i1 + 1], d2, color=ORANGE, lw=0.9, alpha=0.7)
        ax2.set_ylabel("Curvature d$^2$z/dx$^2$ (m$^{-1}$)", color=ORANGE, fontsize=8)
        ax2.tick_params(axis="y", colors=ORANGE, labelsize=7)
        ax2.spines["top"].set_visible(False)
        ax.fill_between(x[v], z[v], -1, color=FILL, lw=0)
        ax.plot(x[v], z[v], color=PROFILE, lw=1.1, label=f"Profile, transect {ex}")
        ax.axhline(MHW, color=WATER, lw=0.7)
        ax.axhline(BERM_M, color=CREST, lw=0.6, ls=(0, (2, 2)))
        ax.text(VIEW_M[0] + 2, BERM_M + 0.12, f"Berm ({BERM_M:g} m)", fontsize=7, color=CREST)
        ax.plot(0, z[i0], "o", ms=5, color=WATER, zorder=6, label="Shoreline")
        ax.plot(x[i1], z[i1], "^", ms=7, color=CREST, zorder=6, label="Foredune crest")
        k = pick_c[ex]
        if k is not None:
            ax.plot(x[k], z[k], "D", ms=7, color=ORANGE, mec="white", mew=0.6, zorder=7,
                    label=f"Maximum curvature toe ({z[k]:.2f} m)")
        ax.plot(x[p["toe"]], z[p["toe"]], "o", ms=7, mfc="white", mec=TOE, mew=1.2, zorder=7,
                label=f"Line method toe ({z[p['toe']]:.2f} m)")
        ax.set_xlim(*VIEW_M)
        ax.set_ylim(-0.5, np.ceil(np.nanmax(z[v]) + 0.8))
        ax.set_xlabel("Distance landward of the shoreline (m)")
        ax.set_ylabel("Elevation (m NAVD88)")
        _legend(ax)

    def pb(ax):                                                  # the domain's five transects
        for i, col in zip(ids, shades):
            s, z, p = profiles[i][:3]
            if p["shore"] is None:
                continue
            x = s - s[p["shore"]]
            ax.plot(x, z, color=col, lw=0.9, zorder=3)
            k = pick_c[i]
            if k is not None:
                ax.plot(x[k], z[k], "D", ms=6, color=col, mec="white", mew=0.6, zorder=6,
                        label=f"{i}: {z[k]:.2f} m")
        ax.axhline(MHW, color=WATER, lw=0.7)
        ax.axhline(BERM_M, color=CREST, lw=0.6, ls=(0, (2, 2)))
        ax.set_xlim(*VIEW_M)
        ax.set_ylim(-0.5, None)
        ax.set_xlabel("Distance landward of the shoreline (m)")
        ax.set_ylabel("Elevation (m NAVD88)")
        fs.open_frame(ax)
        _legend(ax, title="Transect: curvature toe", title_fontsize=7)

    def pc(ax):                                                  # the domain toe
        r = t[t.gis == gis].set_index("transect")
        xs = np.arange(len(ids))
        ax.scatter(xs, [r.curv_toe_m_navd88.get(i, np.nan) for i in ids], s=34, marker="D", color=ORANGE,
                   zorder=5, label="Curvature toe")
        ax.scatter(xs, [r.toe_m_navd88[i] if r.kept[i] else np.nan for i in ids], s=28, facecolors="white",
                   edgecolors=TOE, zorder=4, label="Line method toe")
        ax.axhline(cmed[gis], color=ORANGE, lw=1.4, label=f"Curvature median ({cmed[gis]:.2f} m)")
        ax.axhline(dom.set_index("gis").toe_m_navd88[gis], color=TOE, lw=1.2, ls=(0, (4, 2)),
                   label=f"Line method median ({dom.set_index('gis').toe_m_navd88[gis]:.2f} m)")
        vals = list(r.curv_toe_m_navd88.dropna()) + list(r.toe_m_navd88[r.kept])
        ax.set_ylim(min(vals) - 0.4, max(vals) + 1.2)
        ax.set_xlim(-0.5, len(ids) - 0.5)
        ax.set_xticks(xs, [str(i) for i in ids])
        ax.set_xlabel(f"Transect (GIS {gis}, south → north)")
        ax.set_ylabel("Dune toe (m NAVD88)")
        fs.open_frame(ax)
        _legend(ax, ncol=2)

    def pd_(ax):                                                 # alongshore
        ax.plot(cmed.index, cmed, color=ORANGE, lw=1.4, zorder=5, label="Maximum curvature (domain median)")
        ax.plot(dom.gis, dom.toe_m_navd88, color=TOE, lw=1.2, ls=(0, (4, 2)), zorder=4,
                label="Line method (domain median)")
        ax.set_xlim(0.5, 90.5)
        lo_, hi_ = min(cmed.min(), dom.toe_m_navd88.min()), max(cmed.max(), dom.toe_m_navd88.max())
        ax.set_ylim(np.floor(lo_ - 0.3), np.ceil(hi_ + 1.0))
        ax.axvline(gis, color=fs.INK_MUTED, lw=0.6, ls=":")
        ax.set_xlabel(getattr(fs, "DOMAIN_AXIS_LABEL", "GIS domain (south → north)"))
        ax.set_ylabel("Dune toe (m NAVD88)")
        fs.open_frame(ax)
        _legend(ax, ncol=2)

    out = M_CURVATURE
    dd = (cmed - dom.set_index("gis").toe_m_navd88).dropna()
    caps = {
        "a": (f"Transect {ex} (GIS {gis}): the maximum-curvature toe (diamond; the largest d2z/dx2 between the "
              f"shoreline and the foredune crest, profile smoothed over {CURV_MAIN_M} m, orange line on the right "
              f"axis, no lower than the berm) and the line method's toe (open circle)."),
        "b": f"The five transects of GIS {gis}, aligned at their shorelines, with their maximum-curvature toes.",
        "c": (f"The GIS {gis} domain toe: the median of the curvature toes (solid) against the line method's "
              f"(dashed); open circles, the line method's transect toes."),
        "d": (f"Domain toes alongshore, maximum curvature against the line method (median of each domain's "
              f"transects). Curvature minus line: mean {dd.mean():+.2f} m, RMSE {np.sqrt((dd ** 2).mean()):.2f} m, "
              f"more than 0.5 m apart in {int((dd.abs() > 0.5).sum())} of {len(dd)} domains."),
    }
    for key, draw in zip("abcd", (pa, pb, pc, pd_)):
        fig, ax = plt.subplots(figsize=fs.figsize("double", height=2.9), layout="constrained")
        draw(ax)
        png = out / f"stockdon_{key}.png"
        fs.save(fig, png, close=True)
        fs.record_caption(png, "Maximum-curvature dune toe (Stockdon et al. 2007), 1996 lidar. " + caps[key])
        print(f"wrote {png}")
    fig = plt.figure(figsize=(fs.figsize("double")[0], 10.0), layout="constrained")
    gs = fig.add_gridspec(4, 1)
    for k, (key, draw) in enumerate(zip("abcd", (pa, pb, pc, pd_))):
        ax = fig.add_subplot(gs[k])
        draw(ax)
        ax.text(-0.02, 1.02, key, transform=ax.transAxes, fontsize=10, fontweight="bold", ha="right", va="bottom")
    png = out / f"stockdon_method_{YEAR}.png"
    fs.save(fig, png, close=True)
    fs.record_caption(png, "The dune toe by maximum curvature (Stockdon et al. 2007), 1996 lidar. " + " ".join(
        f"({k}) {v}" for k, v in caps.items()))
    print(f"wrote {png}")


def curvature_example(t, dom, profiles, *_):
    """Why the maximum-curvature toe fails here, simply: (left) an idealised profile, a flat beach meeting a steep
    dune in one sharp bend, where the largest d2z/dx2 is the toe; (right) transect 119 at GIS 22, a beach that rises
    gradually into the dune with small kinks on the face, where the largest d2z/dx2 is a kink halfway up. Top: the
    profile and both toes; bottom: the curvature. Writes figures/curvature_why.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer import hat_figure_style as fs
    fs.apply_style()
    ORANGE = "#E07B39"
    # the idealised profile: 1:30 beach to a 2 m toe at 30 m, a 1:3 dune face to a 6 m crest, a rounded top
    xi = np.arange(0, 61, 1.0)
    zi = np.where(xi <= 30, MHW + (2.0 - MHW) * xi / 30, np.minimum(2.0 + (xi - 30) / 3, 6.0))
    zi = pd.Series(zi).rolling(3, center=True, min_periods=1).mean().to_numpy()
    ids = list(t[t.gis == EXAMPLE_GIS].transect)
    ex = ids[len(ids) // 2]
    s, z, p = profiles[ex][:3]
    i0, i1 = p["shore"], p["crest"]
    xr, zr = s[i0:i1 + 1] - s[i0], z[i0:i1 + 1]

    def toes(x, zz):
        sm = pd.Series(zz).rolling(CURV_MAIN_M, center=True, min_periods=1).mean().to_numpy()
        d2 = np.gradient(np.gradient(sm, STEP_M), STEP_M)
        d2c = d2.copy()
        d2c[:2], d2c[-2:] = -np.inf, -np.inf
        kc = int(np.argmax(d2c))
        sag = np.linspace(zz[0], zz[-1], len(zz)) - zz
        kl = int(np.argmax(sag))
        return d2, kc, kl

    fig, axes = plt.subplots(2, 2, figsize=(fs.figsize("double")[0], 4.6), sharex="col",
                             gridspec_kw=dict(height_ratios=(2, 1)), layout="constrained")
    for col, (x, zz, title) in enumerate(((xi, zi, "Idealised: flat beach, steep dune"),
                                          (xr, zr, f"Hatteras: transect {ex}, GIS {EXAMPLE_GIS}, 1996"))):
        d2, kc, kl = toes(x, zz)
        ax, axc = axes[0, col], axes[1, col]
        ax.fill_between(x, zz, min(zz.min(), 0) - 0.5, color=FILL, lw=0)
        ax.plot(x, zz, color=PROFILE, lw=1.2)
        ax.plot([x[0], x[-1]], [zz[0], zz[-1]], color=CHORD, lw=0.7, ls=(0, (4, 2)))
        ax.plot(x[kl], zz[kl], "o", ms=8, mfc="white", mec=TOE, mew=1.4, zorder=6,
                label=f"Line method ({zz[kl]:.1f} m)")
        ax.plot(x[kc], zz[kc], "D", ms=6, color=ORANGE, mec="white", mew=0.6, zorder=7,
                label=f"Maximum curvature ({zz[kc]:.1f} m)")
        ax.set_title(title, loc="left", fontsize=9)
        ax.set_ylabel("Elevation (m NAVD88)" if col == 0 else "")
        ax.set_ylim(min(zz.min(), 0) - 0.5, zz.max() + 2.2)
        fs.open_frame(ax)
        _legend(ax)
        axc.axhline(0, color=fs.INK_MUTED, lw=0.5)
        axc.plot(x, d2, color=ORANGE, lw=1.0)
        axc.plot(x[kc], d2[kc], "D", ms=5, color=ORANGE, mec="white", mew=0.5)
        axc.set_ylabel("d$^2$z/dx$^2$ (m$^{-1}$)" if col == 0 else "")
        axc.set_xlabel("Distance landward of the shoreline (m)")
        fs.open_frame(axc)
        for a, k in ((ax, "ab"[col]), (axc, "cd"[col])):
            a.text(-0.02, 1.02, k, transform=a.transAxes, fontsize=10, fontweight="bold", ha="right", va="bottom")
    axes[0, 1].annotate("sharpest bend:\na kink on the dune face", xy=(xr[toes(xr, zr)[1]], zr[toes(xr, zr)[1]]),
                        xytext=(-95, 18), textcoords="offset points", fontsize=7, color=ORANGE,
                        arrowprops=dict(arrowstyle="->", color=ORANGE, lw=0.7))
    axes[0, 0].annotate("one sharp bend:\nthe toe", xy=(xi[toes(xi, zi)[1]], zi[toes(xi, zi)[1]]), xytext=(-75, 22),
                        textcoords="offset points", fontsize=7, color=ORANGE,
                        arrowprops=dict(arrowstyle="->", color=ORANGE, lw=0.7))
    png = M_CURVATURE / "curvature_why.png"
    fs.save(fig, png, close=True)
    fs.record_caption(png, (
        f"Why the maximum-curvature toe (Stockdon et al. 2007) fails on these profiles. (a, c) An idealised profile, "
        f"a 1:30 beach meeting a 1:3 dune face in one sharp bend: the largest curvature (d2z/dx2, the profile "
        f"smoothed over {CURV_MAIN_M} m) is the toe, where the line method also puts it. (b, d) Transect {ex} "
        f"(GIS {EXAMPLE_GIS}, 1996 lidar), shoreline to foredune crest: the beach rises gradually into the dune, so "
        f"the toe is a broad, gentle bend, and the largest curvature is a small kink halfway up the dune face. The "
        f"line method (open circle; the point furthest below the shoreline-to-crest line, dashed) still finds the "
        f"foot of the slope. Berm floor not applied here."))
    print(f"wrote {png}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gis", type=int, default=EXAMPLE_GIS, help="the example domain")
    ap.add_argument("--check", action="store_true", help="the check against the line method")
    ap.add_argument("--why", action="store_true", help="why maximum curvature fails, one example")
    a = ap.parse_args()
    if a.check:
        return curvature_check(*run())
    if a.why:
        return curvature_example(*run())
    stockdon_figures(*run(), gis=a.gis)


if __name__ == "__main__":
    main()
