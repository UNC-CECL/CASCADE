"""
Checks on our dune toe (2026-11-05): which crest (the front dune inside the model's dune search window, the window's
highest point, or the front dune within 120 m), and how much the foredune dip that defines the front dune matters
(0.5, 1 and 1.5 m).

    python scripts/real-world-cases/Hatteras/input_prep/9-dune-toes/HAT_toe_checks.py --crest-choice [--gis 22]
    python scripts/real-world-cases/Hatteras/input_prep/9-dune-toes/HAT_toe_checks.py --dip

Writes, in data/hatteras_init/9-dune-toes/1-transects-1996/4-checks/ (D:): crest-choice/ (an example transect and
the domain toes with each crest) and dip/ (dip_sensitivity_1996.csv and its figure).

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
from HAT_toe_common import (BERM_M, CHECKS, CREST_SEARCH_M, EXAMPLE_GIS, FOREDUNE_DROP_M, MHW, MIN_KEPT,
    SHORE_RUN_M, TOE, TOE_BELOW_CREST_M, TOE_MIN_M_NAVD88, WATER, WINDOW, YEAR, _legend)  # noqa: E402
from HAT_toe_method_line import (_Ctx, _medians, _profile, _row, pick, run)  # noqa: E402


def crest_choice(t, dom, profiles, tr, domains, src, scale, offset, gis):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer import hat_figure_style as fs
    fs.apply_style()
    c = _Ctx(t, dom, profiles, tr, domains, src, scale, offset, gis)
    dd = dom.toe_m_navd88 - dom.current_toe_m_navd88
    fmt = dict(gis=gis, ex=c.ex, offset=offset, mhw=MHW, run=SHORE_RUN_M, win=CREST_SEARCH_M, berm=BERM_M,
               drop=FOREDUNE_DROP_M, tmin=TOE_MIN_M_NAVD88, below=TOE_BELOW_CREST_M, nmin=MIN_KEPT,
               nmed=int((dom.source == "median").sum()), nfill=int((dom.source == "filled").sum()),
               dmean=dd.mean(), drmse=float(np.sqrt((dd ** 2).mean())))

    # The two crest definitions compared: the front dune (first peak with a 1 m dip) and the model's window
    dw = dom.toe_window_m_navd88 - dom.current_toe_m_navd88
    dfw = dom.toe_window_m_navd88 - dom.toe_m_navd88
    fmt |= dict(wmean=dw.mean(), wrmse=float(np.sqrt((dw ** 2).mean())), fwmean=dfw.mean(),
                nbig=int((dfw.abs() > 0.5).sum()), nbehind=int((t.window_behind_front_dune_m > 10).sum()),
                nt=len(t), nw=int(t.window_kept.sum()), nf=int(t.kept.sum()))
    cmp_caps = (
        ("compare_example", panel_compare_example,
         "Transect {ex} (GIS {gis}) with the two crest definitions: the front dune (blue; the most seaward peak "
         "with a {drop:g} m dip behind it, within {win:.0f} m of the shoreline) and the highest point of the "
         "model's dune search window (teal band; the window the model's dune crest is taken from), each with the "
         "toe under its shoreline-to-crest chord."),
        ("compare_alongshore", panel_compare_alongshore,
         "Domain dune toes from the 1996 transects with the two crest definitions (median of each domain's retained "
         "transects: front dune {nf}/{nt} retained, model window {nw}/{nt}). Window minus front dune: mean "
         "{fwmean:+.2f} m, more than 0.5 m apart in {nbig} domains; in {nbehind} transects the window's highest "
         "point is a ridge more than 10 m behind the front dune."),
    )
    for name, draw, cap in cmp_caps:
        fig, ax = plt.subplots(figsize=fs.figsize("double", height=2.9), layout="constrained")
        draw(ax, c, fs, fig)
        png = CHECKS / "crest-choice" / f"{name}.png"
        fs.save(fig, png, close=True)
        fs.record_caption(png, cap.format(**fmt))
        print(f"wrote {png}")
    fig, axes = plt.subplots(2, 1, figsize=(fs.figsize("double")[0], 5.6), layout="constrained")
    for k, (ax, (name, draw, cap)) in enumerate(zip(axes, cmp_caps)):
        draw(ax, c, fs, fig)
        ax.text(-0.02, 1.02, "ab"[k], transform=ax.transAxes, fontsize=10, fontweight="bold", ha="right",
                va="bottom")
    png = CHECKS / "crest-choice" / f"compare_crest_methods_{YEAR}.png"
    fs.save(fig, png, close=True)
    fs.record_caption(png, "Two crest definitions for the dune toe. " + " ".join(
        f"({'ab'[k]}) " + cap.format(**fmt) for k, (_, _, cap) in enumerate(cmp_caps)))
    print(f"wrote {png}")


def panel_compare_example(ax, c, fs, fig=None):
    _profile(ax, c, fs)
    s, z, p, _, pw, wa = c.profiles[c.ex][:6]
    o = s[p["shore"]]
    if wa:
        ax.axvspan(wa[1] - o, wa[2] - o, color=WINDOW, alpha=0.12, lw=0, zorder=0, label="Model dune window")
    ax.plot(0, z[p["shore"]], "o", ms=5, color=WATER, zorder=6, label="Shoreline")
    for q, col, name, dash in ((p, TOE, "Front dune in window (used)", (0, (4, 2))),
                               (pw, WINDOW, "Window highest point", (0, (1.5, 1.5)))):
        if q["crest"] is None:
            continue
        cx, zc = s[q["crest"]] - o, z[q["crest"]]
        ax.plot([0, cx], [z[p["shore"]], zc], color=col, lw=0.8, ls=dash, zorder=3)
        ax.plot(cx, zc, "^", ms=7, color=col, zorder=6, label=f"{name}: crest {zc:.2f} m")
        if q["toe"] is not None:
            ax.plot(s[q["toe"]] - o, z[q["toe"]], "o", ms=7, mfc=col, mec="white", mew=0.8, zorder=7,
                    label=f"{name}: toe {z[q['toe']]:.2f} m")
    _legend(ax)


def panel_compare_alongshore(ax, c, fs, fig=None):
    dom = c.dom
    ax.plot(dom.gis, dom.toe_free_m_navd88, color="#9ECAE1", lw=1.2, ls=(0, (4, 2)), zorder=3,
            label="Front dune within 120 m")
    ax.plot(dom.gis, dom.toe_window_m_navd88, color=WINDOW, lw=1.4, zorder=4, label="Window highest point")
    ax.plot(dom.gis, dom.toe_m_navd88, color=TOE, lw=1.6, zorder=5, label="Front dune in window (used)")
    ax.set_xlim(0.5, 90.5)
    v = pd.concat([dom.toe_m_navd88, dom.toe_window_m_navd88, dom.toe_free_m_navd88])
    ax.set_ylim(np.floor(v.min() - 0.3), np.ceil(v.max() + 1.0))
    ax.set_xlabel(getattr(fs, "DOMAIN_AXIS_LABEL", "GIS domain (south → north)"))
    ax.set_ylabel("Dune toe (m NAVD88)")
    fs.open_frame(ax)
    _legend(ax, ncol=3)


DIP_TEST_M = (0.5, 1.0, 1.5)


def dip_test(t, dom, profiles, *_):
    """The foredune dip threshold: the same profiles re-picked at each DIP_TEST_M; domain medians and how many
    domains move against the 1 m default. Writes dip_sensitivity_<year>.csv and figures/dip_sensitivity.png."""
    rows = []
    for drop in DIP_TEST_M:
        for r in t.itertuples():
            s, z = profiles[r.transect][:2]
            wa = profiles[r.transect][5]
            p = pick(s, z, wa[1:] if wa else None, drop=drop)
            rows.append(dict(drop_m=drop, transect=r.transect, gis=r.gis,
                             **_row(s, z, p)))
    d = pd.DataFrame(rows)
    out = []
    for drop, g in d.groupby("drop_m"):
        med, n, src = _medians(g)
        out.append(pd.DataFrame({"drop_m": drop, "gis": range(1, 91), "toe_m_navd88": med, "n_kept": n,
                                 "source": src}))
    m = pd.concat(out, ignore_index=True)
    w = m.pivot(index="gis", columns="drop_m", values="toe_m_navd88")
    (CHECKS / "dip").mkdir(parents=True, exist_ok=True)
    w.round(3).add_prefix("toe_dip_").to_csv(CHECKS / "dip" / f"dip_sensitivity_{YEAR}.csv")
    base = w[1.0]
    for drop in DIP_TEST_M:
        g = d[d.drop_m == drop]
        diff = w[drop] - base
        moved = (g.set_index("transect").crest_s_m - d[d.drop_m == 1.0].set_index("transect").crest_s_m).abs() > 1
        print(f"dip {drop:g} m: kept {int(g.kept.sum())}/{len(g)}, crest moved on {int(moved.sum())} transects; "
              f"domain toe minus 1 m: mean {diff.mean():+.2f} m, |diff| > 0.25 m in {int((diff.abs() > 0.25).sum())}"
              f", > 0.5 m in {int((diff.abs() > 0.5).sum())} domains"
              + (f" ({', '.join(str(i) for i in diff.index[diff.abs() > 0.5])})" if (diff.abs() > 0.5).any() else ""))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer import hat_figure_style as fs
    fs.apply_style()
    fig, ax = plt.subplots(figsize=fs.figsize("double", height=2.9), layout="constrained")
    for drop, col, lw in zip(DIP_TEST_M, ("#9ECAE1", "#08519C", "#4292C6"), (1.2, 1.6, 1.2)):
        ax.plot(w.index, w[drop], color=col, lw=lw, ls="-" if drop == 1.0 else (0, (4, 2)),
                zorder=5 if drop == 1.0 else 4, label=f"{drop:g} m dip" + (" (default)" if drop == 1.0 else ""))
    ax.set_xlim(0.5, 90.5)
    ax.set_xlabel(getattr(fs, "DOMAIN_AXIS_LABEL", "GIS domain (south → north)"))
    ax.set_ylabel("Dune toe (m NAVD88)")
    ax.set_ylim(np.floor(w.min().min() - 0.3), np.ceil(w.max().max() + 1.0))
    fs.open_frame(ax)
    _legend(ax, ncol=3)
    png = CHECKS / "dip" / "dip_sensitivity.png"
    fs.save(fig, png, close=True)
    fs.record_caption(png, (
        f"Sensitivity of the domain dune toe (median of the retained transects, {YEAR} lidar) to the dip that "
        f"marks the foredune: the most seaward peak with at least 0.5, 1.0 (default) or 1.5 m of dip behind it, "
        f"within {CREST_SEARCH_M:.0f} m of the shoreline."))
    print(f"wrote {png}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gis", type=int, default=EXAMPLE_GIS, help="the example domain")
    ap.add_argument("--crest-choice", action="store_true", help="the three crest definitions compared")
    ap.add_argument("--dip", action="store_true", help="the foredune dip threshold test")
    a = ap.parse_args()
    if a.dip:
        dip_test(*run())
    if a.crest_choice or not a.dip:
        crest_choice(*run(), a.gis)


if __name__ == "__main__":
    main()
