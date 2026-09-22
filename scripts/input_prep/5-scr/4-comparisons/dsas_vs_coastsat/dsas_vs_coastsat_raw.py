"""
dsas_vs_coastsat_raw.py
==============================================================================
The two shoreline-rate sources against each other with NO SMOOTHING: the raw
per-domain mean LRR from DSAS and from CoastSat, on the two DSAS windows
(Hannah, 2026-09-22: "I want to see it without smoothing").

    4-comparisons/dsas_vs_coastsat/calendar_windows/dsas_vs_coastsat_raw.png
    4-comparisons/dsas_vs_coastsat/calendar_windows/slides/ (--slide)
    4-comparisons/dsas_vs_coastsat/calendar_windows/supporting/*.csv

WHY IT IS SEPARATE FROM 6-scr-smooth/dsas_vs_coastsat/
    That folder exists to argue about the SMOOTHING -- every figure in it
    draws a LOESS, and its own README calls these windows retired. The
    question here is different and prior to it: before any smoothing, do the
    two sources say the same thing about the same 500 m of beach? So it sits
    with the other comparisons, and the figure carries no LOESS at all.

WHAT THE CoastSat SIDE IS
    REFIT HERE from the current time series (Hannah, 2026-09-22), not read
    from the archive: every transect in the current transect_domain_lookup.csv
    is fitted over the DSAS window with coastsat_lrr.compute_lrr --
    the same loader, the same date filter, the same OLS, the same 3-position
    minimum as the live windows -- and averaged per domain. Only the window
    differs from 3-rates/coastsat/lrr/.

    The refit is NOT written to 3-rates/coastsat/lrr/. hat_observed_rates
    .windows() enumerates that tree by scanning it, so a 1978_1997 folder
    there would become a window for every caller: rates_figures would draw
    figures for it, and the shared y bound of EVERY window figure is the
    largest domain mean over all windows, so the existing figures would
    change. These years are a comparison, not a rate product, and they stay
    inside this folder.

    The archived fits from 5-scr/archive/coastsat_lrr_superseded_20260810/
    are still read, for one purpose: the table carries them beside the refit
    so the two can be differenced. Nothing in this folder feeds a run.

    READ THE FIRST PANEL WITH CARE. CoastSat imagery begins in 1984, so its
    "1978-1997" rate is fitted from 1984-06-17 -- 13.5 years against DSAS's
    19, and it misses the six years at the start entirely. The second window
    is a fair comparison (CoastSat 1997-01-12 to 2019-12-28); the first is
    two different periods with one label, which is worth more of the
    disagreement than any method difference.

SLIDE VERSION
    --slide draws the same two panels on a 3.4 in canvas, as
    dsas_vs_coastsat_raw_slide.png. The reach is NOT cut down the way
    lrr_transect_zoom's slide is: the whole island IS the comparison here, and
    the thing a viewer reads off it -- the two lines apart in (a), together in
    (b) -- survives the smaller canvas, where a nine-domain crop would not
    show it at all. Only the labels, ticks and line weight come down.

USAGE
    python scripts/input_prep/5-scr/4-comparisons/dsas_vs_coastsat/dsas_vs_coastsat_raw.py
    python scripts/input_prep/5-scr/4-comparisons/dsas_vs_coastsat/dsas_vs_coastsat_raw.py --slide
==============================================================================
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_REPO = next(_p for _p in Path(__file__).resolve().parents
             if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(_REPO / "scripts"))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import MultipleLocator  # noqa: E402

sys.path.insert(0, str(_REPO / "scripts" / "input_prep" / "5-scr" / "lib"))
import scr_paths  # noqa: E402,F401  (5-scr sibling modules onto sys.path)
from coastsat_lrr import (  # noqa: E402
    _empty_lrr, compute_lrr, filter_dates, load_timeseries,
)
from site_layer import hat_observed_rates as obs  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C_1984, C_1997, DOMAIN_AXIS_LABEL, INK_MUTED, apply_style, caption,
    figsize, open_frame, save, structures, support_dir, town_bands,
)

N = 90
WINDOWS = ((1978, 1997), (1997, 2019))
OUT = obs.COMPARISONS / "dsas_vs_coastsat" / "calendar_windows"
STEM = "dsas_vs_coastsat_raw"
C_DSAS, C_CS = C_1984, C_1997
LW = 1.1
SLIDE_W_IN = 3.4
SLIDE_H_IN = 4.2
LW_SLIDE = 0.8
MIN_OBS = 3   # as the live windows, coastsat_domain_lrr_fixed.MIN_OBS


def dsas(start, end):
    """Per-domain mean LRR from the DSAS rate table, on domains 1..90."""
    df = pd.read_csv(obs.DSAS_ROOT / f"dsas_{start}_{end}_rates.csv")
    df = df.rename(columns={"domain_id": "domain", "MEAN_LRR": "lrr"})
    return (df[["domain", "lrr"]].groupby("domain")["lrr"].mean()
            .reindex(range(1, N + 1)))


def coastsat_archived(start, end):
    """Per-domain mean LRR as the retired fit recorded it, for the diff only."""
    df = pd.read_csv(obs.COASTSAT_LRR_SUPERSEDED / f"{start}_{end}"
                     / "domain_lrr_summary.csv")
    df = df.rename(columns={"domain_number": "domain", "mean_lrr": "lrr"})
    return (df[["domain", "lrr"]].groupby("domain")["lrr"].mean()
            .reindex(range(1, N + 1)))


def _series_cache():
    """{transect_id: chainage frame} for every transect in the current lookup.

    Read once and fitted to both windows, rather than once per window: the
    time series are the same file either way.
    """
    lookup = pd.read_csv(obs.transect_lookup())
    lookup = lookup[lookup["domain_number"].between(1, N)]
    cache, missing = {}, 0
    for tid in lookup["transect_id"]:
        site = tid.rsplit("_", 1)[0]
        path = obs.COASTSAT_TIMESERIES / f"{site}_timeseries" / f"{tid}.csv"
        if not path.is_file():
            missing += 1
            continue
        cache[tid] = load_timeseries(str(path))
    if missing:
        print(f"  {missing} transects in the lookup have no time series")
    return lookup, cache


def coastsat_refit(start, end, lookup, cache):
    """Per-domain mean LRR refitted from the current time series.

    Same method as the live windows: every position inside the calendar
    window, no outlier filter, no weighting, at least MIN_OBS positions.
    Returns (per-domain Series, per-transect frame).
    """
    rows = []
    for tid, dom in zip(lookup["transect_id"], lookup["domain_number"]):
        df = cache.get(tid)
        if df is None:
            continue
        sel = filter_dates(df, f"{start}-01-01", f"{end}-12-31")
        r = compute_lrr(sel) if len(sel) >= MIN_OBS else _empty_lrr(len(sel))
        rows.append(dict(transect_id=tid, domain_number=int(dom), **r))
    t = pd.DataFrame(rows)
    per_domain = (t.dropna(subset=["lrr_m_yr"])
                  .groupby("domain_number")["lrr_m_yr"].mean()
                  .reindex(range(1, N + 1)))
    return per_domain, t


def agreement(a, b):
    """n, bias, RMSE and r over the domains where both sources have a value.

    bias is CoastSat minus DSAS, so positive means CoastSat reports the more
    seaward rate.
    """
    ok = a.notna() & b.notna()
    d = (b[ok] - a[ok]).to_numpy(float)
    r = (float(np.corrcoef(a[ok], b[ok])[0, 1]) if ok.sum() > 2 else float("nan"))
    return dict(n=int(ok.sum()), bias=float(d.mean()),
                rmse=float(np.sqrt((d ** 2).mean())), r=r,
                max_abs=float(np.abs(d).max()),
                max_at=int((b - a).abs().idxmax()))


def figure(series, stats, slide=False):
    x = np.arange(1, N + 1)
    lw = LW_SLIDE if slide else LW
    size = (figsize(SLIDE_W_IN, height=SLIDE_H_IN) if slide
            else figsize("double", height=5.0))
    fig, axes = plt.subplots(len(WINDOWS), 1, sharex=True, sharey=True,
                             figsize=size, constrained_layout=True)
    for i, ((s, e), ax) in enumerate(zip(WINDOWS, axes)):
        ax.set_xlim(0.5, N + 0.5)
        # Village names and structure labels are page furniture: at 3.4 in
        # they collide with each other and with the data (2026-09-22). The
        # bands stay, unlabelled, so the villages are still locatable.
        town_bands(ax, label=(i == 0 and not slide))
        ax.axhline(0, color=INK_MUTED, lw=0.6, zorder=2)
        ax.plot(x, series[(s, e)]["dsas"], color=C_DSAS, lw=lw, zorder=4)
        ax.plot(x, series[(s, e)]["cs"], color=C_CS, lw=lw, zorder=5)
        structures(ax, i == 0 and not slide, 4.0)
        ax.yaxis.set_major_locator(MultipleLocator(2))
        ax.xaxis.set_major_locator(MultipleLocator(20 if slide else 10))
        ax.xaxis.set_minor_locator(MultipleLocator(10 if slide else 5))
        ax.yaxis.grid(True, zorder=0)
        ax.set_axisbelow(True)
        ax.set_title(f"({chr(97 + i)}) {s}–{e}", loc="left", fontsize=8)
        open_frame(ax)
    axes[-1].set_xlabel("GIS domain (S → N)" if slide else DOMAIN_AXIS_LABEL)
    fig.supylabel("Shoreline change rate (m/yr)", fontsize=9)
    fig.legend([Line2D([], [], color=C_DSAS, lw=lw),
                Line2D([], [], color=C_CS, lw=lw)],
               (["DSAS", "CoastSat"] if slide else
                ["DSAS (digitized shorelines)", "CoastSat (satellite)"]),
               loc="outside lower center", ncol=2, frameon=False)

    t = "; ".join(
        f"{s}–{e} bias {stats[(s, e)]['bias']:+.2f}, RMSE "
        f"{stats[(s, e)]['rmse']:.2f} m/yr, r {stats[(s, e)]['r']:.2f} "
        f"(n {stats[(s, e)]['n']})" for s, e in WINDOWS)
    caption(fig, (
        "Observed shoreline change rate by GIS domain (1 at Cape Point, 90 at "
        "Pea Island) from the two sources, UNSMOOTHED: each line is the plain "
        "per-domain mean rate, DSAS from the digitized shoreline transects and "
        "CoastSat from the satellite transects, seaward positive. No LOESS is "
        "applied to either — this is the comparison before any smoothing, so "
        "the per-domain disagreement is shown at full amplitude. Agreement, "
        f"CoastSat minus DSAS: {t}. "
        "The first panel is not a like-for-like window: CoastSat imagery "
        "begins in 1984, so its 1978–1997 rate is fitted from 1984-06-17 over "
        "13.5 years against the DSAS 19, missing the first six entirely. The "
        "second panel compares the same years on both sides (CoastSat "
        "1997-01-12 to 2019-12-28). The CoastSat rates are refitted here from "
        "the current time series and transect lookup, by the same method as "
        "the live windows, because no window in 3-rates/coastsat/lrr/ covers "
        "these years; they are a comparison, not a scoring target, and are "
        "not written into that tree. Village spans are shaded; the solid hairline is the "
        "Buxton groin and the dotted hairlines are the Avon and Rodanthe "
        "piers. Both panels share one y axis."))
    return fig


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--slide", action="store_true",
                    help=f"{SLIDE_W_IN:g} in canvas, to <stem>_slide.png")
    slide = ap.parse_args(argv).slide
    apply_style()
    lookup, cache = _series_cache()
    print(f"  {len(cache)} time series loaded")
    series, stats, cols = {}, {}, {"domain": np.arange(1, N + 1)}
    for s, e in WINDOWS:
        a, _t = dsas(s, e), None
        b, t = coastsat_refit(s, e, lookup, cache)
        old = coastsat_archived(s, e)
        series[(s, e)] = {"dsas": a, "cs": b}
        stats[(s, e)] = agreement(a, b)
        gap = (b - old).abs()
        print(f"  {s}-{e}: refit vs archived fit, mean |diff| {gap.mean():.3f}"
              f"  max {gap.max():.3f} m/yr  (n transects fitted "
              f"{int(t['lrr_m_yr'].notna().sum())})")
        cols[f"dsas_{s}_{e}"] = a.to_numpy(float)
        cols[f"coastsat_{s}_{e}"] = b.to_numpy(float)
        cols[f"difference_{s}_{e}"] = (b - a).to_numpy(float)
        cols[f"coastsat_archived_{s}_{e}"] = old.to_numpy(float)

    fig = figure(series, stats, slide=slide)
    folder = OUT / "slides" if slide else OUT
    out = save(fig, folder / (STEM + ("_slide" if slide else "")))
    plt.close(fig)
    table = support_dir(OUT) / f"{STEM}.csv"
    pd.DataFrame(cols).round(3).to_csv(table, index=False)

    for (s, e), v in stats.items():
        print(f"{s}-{e}: n {v['n']}  bias {v['bias']:+.2f}  RMSE {v['rmse']:.2f}"
              f"  r {v['r']:.2f}  largest gap {v['max_abs']:.2f} m/yr at GIS {v['max_at']}")
    for p in out + [table]:
        print(f"wrote    {p.relative_to(_REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
