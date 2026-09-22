"""
coastsat_lrr_transect_zoom.py
==============================================================================
One window's LRR at TRANSECT resolution over a short reach, with the transect
itself on the x axis (Hannah, 2026-09-22: "x being the transects instead of
the domains").

    coastsat/lrr/<w>/transects/lrr_transects_<w>_gis<lo>-<hi>.png
    coastsat/lrr/<w>/transects/slides/  the --slide versions

    Its own subfolder because the family grows a file per reach and per
    variant, and the window folder's convention is tables plus the one figure
    for the window (2026-09-22). smoothing_windows_<w>.png stays at the top
    level with lrr_<w>.png: coastsat_total_change.py and five PROVENANCE.md
    files cross-reference it by that path.

WHY IT EXISTS
    Every other figure in 3-rates puts the GIS domain on x and shows the
    transects as a cloud behind the domain mean. That is the right axis for
    comparing against the model, whose cell IS the domain, but it hides the
    step the averaging actually takes: which transects fall in which domain,
    how far apart they sit, and how much of the domain mean is one transect.
    This figure is that step, drawn -- the example for explaining the chain,
    not a product any run reads.

WHAT IT SHOWS
    One marker per CoastSat transect, in order south to north, coloured by its
    own sign (the house blue / red pair) with its 95 % regression confidence
    interval as a whisker. The domain each transect belongs to is the shaded
    band behind it, labelled along the top; the flat segment across each band
    is that domain's mean -- the plain average of the markers above it, the
    value domain_lrr_summary.csv holds. That is the whole figure by default.

    --target adds the scoring curve as a second set of segments, written to
    <stem>_with_target.png so it never overwrites the plain version. It is not
    an average of anything inside the band: the LOESS reads 5 km either way,
    so it can sit off the mean, which is the point of drawing it -- but it is
    a third quantity on the panel and needs explaining before it can be read,
    so the teaching version leaves it off.

THE REACH
    Default GIS 64-72, the Tri-Village gradient: the domain mean crosses zero
    at 64 and reaches +2.8 m/yr at 67, so a real alongshore signal is resolved
    transect by transect rather than a flat stretch where every marker is the
    same number. --gis takes any span.

USAGE
    python scripts/input_prep/5-scr/3-rates/coastsat/lrr/coastsat_lrr_transect_zoom.py
    python scripts/input_prep/5-scr/3-rates/coastsat/lrr/coastsat_lrr_transect_zoom.py --window 1996_2024 --gis 28 36
    python scripts/input_prep/5-scr/3-rates/coastsat/lrr/coastsat_lrr_transect_zoom.py --target
    python scripts/input_prep/5-scr/3-rates/coastsat/lrr/coastsat_lrr_transect_zoom.py --slide
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
sys.path.insert(0, str(_REPO / "scripts" / "input_prep" / "5-scr" / "lib"))
import scr_paths  # noqa: E402,F401  (5-scr sibling modules onto sys.path)

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.legend_handler import HandlerTuple  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import MultipleLocator  # noqa: E402

import coastsat_lrr_windows as cw  # noqa: E402
import rates_figures as rf  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    GRID_C, INK_MUTED, SMOOTH_RAMP, apply_style, caption, figsize, save,
)
from cascade_pipeline.coastsat_loess import (  # noqa: E402
    DEFAULT_LOESS, loess_transect_values,
)
from site_layer.hat_observed_rates import COASTSAT_LRR_ROOT  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM  # noqa: E402

DEFAULT_WINDOW = "1996_2024"
DEFAULT_GIS = (64, 72)
BAND_C = GRID_C          # the alternating domain bands
MEAN_LW = 1.6
TARGET_LW = 1.6
Y_PAD = 0.4              # m/yr of headroom over the whiskers
# --slide: the same figure on a canvas that fits a slide beside text. The
# house font sizes are absolute, so a smaller canvas makes the type
# relatively LARGER, which is what a projected figure needs; only the marker
# size, the tick spacing and the legend columns have to come down with it
# (Hannah, 2026-09-22 -- same layout, just smaller).
SLIDE_W_IN = 3.4
SLIDE_ASPECT = 0.80
# A slide panel is a third the width of the page figure, so the same nine
# domains would put 90 markers in 3.4 in and the gradient would read as a
# smear. --slide therefore narrows the REACH as well as the canvas, to four
# domains across the peak, unless --gis says otherwise (Hannah, 2026-09-22).
SLIDE_GIS = (66, 69)
DOT_S = 11.0
DOT_S_SLIDE = 9.0
XTICK = 10
XTICK_SLIDE = 10
LOESS_LW = 1.6
LOESS_WINDOW = max(DEFAULT_LOESS.window_domains)
MAX_LABELLED_BANDS = 20  # past this, no domain shading and no numbers
# The along-coast axis covers anything from a 4 km reach to the whole 45 km
# island, so neither its unit nor its tick step can be a constant: 1 km ticks
# over the island put 45 labels on top of each other.
MAX_XTICKS = 9
NICE_STEPS_M = (100.0, 200.0, 500.0, 1000.0, 2000.0, 5000.0, 10000.0)
KM_ABOVE_M = 10000.0     # span past which the axis is drawn in km, not m


def load(stem, lo, hi, loess_window=None):
    """Transect rows inside GIS lo..hi, south to north, with an x position.

    x is the transect's ORDER in the reach, not a distance: one unit per
    transect is what puts the transect on the axis. The transects are ~50 m
    apart and evenly spread inside a domain, so the two differ only where a
    domain holds an unusual number of them, and the band edges below are drawn
    from the counts rather than assumed.

    With `loess_window`, a `loess` column carries the smoother's value at each
    transect. It is fitted over the WHOLE island first and sliced to the reach
    afterwards, never fitted to the reach alone: a LOESS reads 5 km either
    way, so a fit stopping at the reach edge would be a different curve from
    the one the target is built through.
    """
    t = pd.read_csv(COASTSAT_LRR_ROOT / stem / "transect_lrr_full.csv")
    t["domain_number"] = t["domain_number"].astype(int)
    t = t[t["domain_number"].between(DOM.first_gis_id, DOM.last_gis_id)].copy()
    t = t.sort_values(["domain_number", "transect_id"]).reset_index(drop=True)
    sp = DOM.domain_spacing_m
    rank = t.groupby("domain_number").cumcount()
    n = t.groupby("domain_number")["domain_number"].transform("count")
    t["along_m"] = ((t["domain_number"] - DOM.first_gis_id) * sp
                    + (rank + 0.5) * (sp / n))
    frac = None
    if loess_window:
        sm, frac = loess_transect_values(t["along_m"].to_numpy(float),
                                         t["lrr_m_yr"].to_numpy(float),
                                         loess_window)
        t["loess"] = sm if sm is not None else np.nan
    t = t[t["domain_number"].between(lo, hi)].reset_index(drop=True)
    t["x"] = np.arange(len(t), dtype=float)
    return t, frac


def along_axis(span_m):
    """(scale, unit, tick step) for an along-coast axis spanning `span_m`.

    scale divides the metre values for plotting, so the numbers on the page
    stay short; the coordinate itself is unchanged.
    """
    scale, unit = ((1000.0, "km") if span_m >= KM_ABOVE_M else (1.0, "m"))
    step = next((v for v in NICE_STEPS_M if span_m / v <= MAX_XTICKS),
                NICE_STEPS_M[-1])
    return scale, unit, step / scale


def bands(t, metres=False):
    """[(domain, left, right)] per domain, in the coordinate on the x axis.

    In metres the edges are the domain's true 500 m boundaries. On the
    transect-order axis there is no such thing -- a domain is however many
    transects fell in it -- so the edges are half a transect either side of
    its first and last.
    """
    out = []
    sp = DOM.domain_spacing_m
    for d in sorted(t["domain_number"].unique()):
        g = t.loc[t["domain_number"] == d]
        if metres:
            left = (d - DOM.first_gis_id) * sp
            right = left + sp
        else:
            left, right = float(g["x"].min()) - 0.5, float(g["x"].max()) + 0.5
        out.append((int(d), float(left), float(right)))
    return out


def figure(stem, lo, hi, with_target=False, slide=False, with_loess=False,
           x="transect", domains_shown=True):
    start, end = (int(v) for v in stem.split("_"))
    metres = (x == "metres")
    t, frac = load(stem, lo, hi, LOESS_WINDOW if with_loess else None)
    if t.empty:
        raise SystemExit(f"no transects in GIS {lo}-{hi} for {stem}")
    xv = t["along_m" if metres else "x"].to_numpy(float)
    xscale, xunit, xstep = along_axis(float(xv.max() - xv.min())) if metres else (1.0, "", 0)
    xv = xv / xscale
    y = t["lrr_m_yr"].to_numpy(float)
    unc = t["unc_m_yr"].to_numpy(float)
    # Off by default (Hannah, 2026-09-22): this figure's job is the transect
    # -> domain step, and the scoring curve is a third thing on the panel
    # that has to be explained before it can be read. --target puts it back,
    # under its own name, so the two versions never overwrite each other.
    target = (rf._target(pd.read_csv(COASTSAT_LRR_ROOT / stem / "transect_lrr_full.csv"))
              if with_target else None)

    size = (figsize(SLIDE_W_IN, aspect=SLIDE_ASPECT) if slide
            else figsize("double", aspect=0.40))
    fig, ax = plt.subplots(figsize=size, constrained_layout=True)

    # Every other domain shaded, so a band is a domain without a legend entry.
    # Past ~20 domains the bands are a comb and the numbers collide, so a long
    # reach gets neither and reads as the plain scatter it is.
    bs = [(d, a / xscale, b / xscale) for d, a, b in bands(t, metres)]
    if domains_shown and len(bs) <= MAX_LABELLED_BANDS:
        for i, (d, x0, x1) in enumerate(bs):
            if i % 2:
                ax.axvspan(x0, x1, color=BAND_C, alpha=0.45, lw=0, zorder=0)
            ax.text((x0 + x1) / 2, 1.0, str(d),
                    transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                    fontsize=6.5, color=INK_MUTED, clip_on=False)

    ax.axhline(0, color=INK_MUTED, lw=0.6, zorder=2)

    col = np.where(y < 0, cw.C_ERODE, cw.C_ACCRETE)
    ax.errorbar(xv, y, yerr=unc, fmt="none", ecolor=INK_MUTED, elinewidth=0.6,
                alpha=0.55, capsize=0, zorder=3)
    ax.scatter(xv, y, s=(DOT_S_SLIDE if slide else DOT_S), c=col,
               linewidths=0, zorder=4)

    # The per-domain values, each flat across the domain it belongs to. The
    # target takes a mid-ramp blue when the LOESS curve is drawn too, so the
    # curve and the domain-resolution version of it are not one colour.
    target_c = SMOOTH_RAMP[1] if (with_loess and target is not None) else SMOOTH_RAMP[-1]
    for d, x0, x1 in (bs if domains_shown else []):
        m = float(t.loc[t["domain_number"] == d, "lrr_m_yr"].mean())
        ax.plot([x0, x1], [m, m], color=INK_MUTED, lw=MEAN_LW,
                solid_capstyle="butt", zorder=5)
        if target is not None and d in target.index and np.isfinite(target.loc[d]):
            ax.plot([x0, x1], [target.loc[d]] * 2, color=target_c, lw=TARGET_LW,
                    solid_capstyle="butt", zorder=6)

    if with_loess and "loess" in t:
        ax.plot(xv, t["loess"].to_numpy(float), color=SMOOTH_RAMP[-1],
                lw=LOESS_LW, solid_capstyle="round", zorder=7)

    if not domains_shown:
        # No bands to align to, so the panel holds the markers and a little
        # air, not the 500 m boundaries they happen to sit between.
        pad = (float(np.median(np.diff(np.sort(xv)))) if len(xv) > 1 else 1.0)
        ax.set_xlim(float(xv.min()) - pad, float(xv.max()) + pad)
    else:
        ax.set_xlim(*((bs[0][1], bs[-1][2]) if metres else (-1.0, len(t))))
    # Bounds over EVERY series drawn, not just the markers: the LOESS reads
    # 5 km beyond the reach and routinely sits below everything in it, so
    # bounds taken from the markers alone would clip the curve.
    drawn = [y - unc, y + unc]
    if domains_shown:
        drawn.append(np.array([float(t.loc[t["domain_number"] == d, "lrr_m_yr"].mean())
                               for d, _, _ in bs]))
    if with_loess and "loess" in t:
        drawn.append(t["loess"].to_numpy(float))
    if target is not None and domains_shown:
        drawn.append(target.reindex([d for d, _, _ in bs]).to_numpy(float))
    ylo = float(np.nanmin(np.concatenate(drawn))) - Y_PAD
    yhi = float(np.nanmax(np.concatenate(drawn))) + Y_PAD
    # The page figure always holds zero -- the line sign is read against is
    # worth the white space in a document. A slide panel is a third the width
    # and cannot spare it, so it crops to the data and the caption says when
    # zero fell off (Hannah, 2026-09-22).
    if not slide:
        ylo, yhi = min(0.0, ylo), max(0.0, yhi)
    ax.set_ylim(ylo, yhi)
    zero_shown = ylo <= 0.0 <= yhi
    if metres:
        ax.xaxis.set_major_locator(MultipleLocator(xstep))
        if len(bs) <= MAX_LABELLED_BANDS:
            ax.xaxis.set_minor_locator(
                MultipleLocator(DOM.domain_spacing_m / xscale))
    else:
        ax.xaxis.set_major_locator(MultipleLocator(XTICK_SLIDE if slide else XTICK))
    ax.yaxis.grid(True, zorder=0)
    ax.set_axisbelow(True)
    banded = domains_shown and len(bs) <= MAX_LABELLED_BANDS
    # The page label does not fit a 3.4 in panel -- it ran off both edges --
    # so the slide keeps only what the axis cannot be read without, and the
    # shading and the domain numbers are left to the caption (2026-09-22).
    if not domains_shown:
        # No domain vocabulary on the panel at all, not even the reach: the
        # figure is the transect fits and nothing else, and where on the
        # island it sits is the caption's job (Hannah, 2026-09-22).
        ax.set_xlabel(f"Along-coast distance ({xunit})" if metres
                      else ("CoastSat transect, S → N" if slide
                            else "CoastSat transect, south → north"))
    elif slide:
        ax.set_xlabel((f"Along-coast distance ({xunit})" if metres
                       else "CoastSat transect, S → N")
                      + f", GIS {lo}–{hi}")
    else:
        ax.set_xlabel(
            (f"Along-coast distance from GIS {DOM.first_gis_id} ({xunit}), "
             "south → north" if metres else "CoastSat transect, south → north")
            + f" (GIS {lo}–{hi})"
            + (", domain shaded and numbered" if banded else ""))
    ax.set_ylabel(cw.Y_LABEL)
    cw.open_frame(ax)

    h = [(Line2D([], [], color=cw.C_ACCRETE, marker="o", ms=3.4, lw=0),
          Line2D([], [], color=cw.C_ERODE, marker="o", ms=3.4, lw=0)),
         Line2D([], [], color=INK_MUTED, lw=0.6)]
    labels = ["Individual transects", "95 % confidence interval"]
    if domains_shown:
        h.append(Line2D([], [], color=INK_MUTED, lw=MEAN_LW))
        labels.append("Domain mean")
    if target is not None:
        h.append(Line2D([], [], color=target_c, lw=TARGET_LW))
        labels.append("Graded target")
    if with_loess:
        h.append(Line2D([], [], color=SMOOTH_RAMP[-1], lw=LOESS_LW))
        labels.append(f"LOESS, {LOESS_WINDOW * 0.5:g} km")
    fig.legend(h, labels, loc="outside lower center",
               ncol=(2 if slide else len(h)), frameon=False,
               handler_map={tuple: HandlerTuple(ndivide=None, pad=0.3)})

    n_dom = len(bs)
    per = t.groupby("domain_number").size()
    caption(fig, (
        f"Observed shoreline change rate at transect resolution, GIS {lo}–{hi} "
        f"({n_dom} domains, {len(t)} CoastSat transects, {per.min()}–{per.max()} "
        f"per domain), {start}–{end}. Each marker is one transect's ordinary-"
        "least-squares slope of shoreline position against date over the calendar "
        f"window (1 January {start} to 31 December {end}), seaward positive, blue "
        "where the shoreline moved seaward and red where it moved landward; the "
        "whisker is that fit's 95 % confidence interval."
        + (" Shaded bands are the 500 m GIS domains, numbered along the top: "
           "the grey segment across each is the plain mean of the markers above "
           "it, the value domain_lrr_summary.csv holds." if domains_shown else
           " Nothing on the panel is aggregated: no domain bands, no domain "
           "means, only the transects as fitted.")
        + (" The dark blue segment is the target a model run is scored against, a "
           "5 km alongshore LOESS of these same transect rates; it is not an "
           "average of the band it sits in — it reads 5 km in both directions, so "
           "it can fall outside the markers beneath it." if target is not None else "")
        + (f" The dark blue curve is the LOESS itself, at the resolution it is "
           f"fitted at: one value per transect, window {LOESS_WINDOW * 0.5:g} km "
           f"(frac {frac:.3f} of the island's transects). It is fitted over the "
           "whole island and sliced to this reach, not fitted to the reach, so "
           "near either edge it is reading transects outside the panel."
           if with_loess else "")
        + (" The x axis is along-coast distance, the coordinate the LOESS "
           f"actually uses, measured from the south end of the modelled reach "
           f"(GIS {DOM.first_gis_id}); transects sit about 50 m apart, so the "
           "spacing on the page is the spacing on the beach."
           if metres else "")
        + " The y axis is set by this reach, not by the island-wide bound of "
        "the window figures"
        + ("" if zero_shown else
           ("; zero is off the panel, so every rate shown is "
            + ("seaward" if ylo > 0 else "landward")))
        + ". "
        + rf._marks_clause(start, end)))
    folder = COASTSAT_LRR_ROOT / stem / "transects" / ("slides" if slide else "")
    name = (f"lrr_transects_{stem}_gis{lo}-{hi}"
            + ("" if domains_shown else "_no_domains")
            + ("_metres" if metres else "")
            + ("_with_loess" if with_loess else "")
            + ("_with_target" if target is not None else "")
            + ("_slide" if slide else ""))
    out = save(fig, folder / name)
    plt.close(fig)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--window", default=DEFAULT_WINDOW, help="e.g. 1996_2024")
    ap.add_argument("--gis", nargs=2, type=int, default=None,
                    metavar=("LO", "HI"),
                    help=f"inclusive GIS domain span (default "
                         f"{DEFAULT_GIS[0]}-{DEFAULT_GIS[1]}, or "
                         f"{SLIDE_GIS[0]}-{SLIDE_GIS[1]} with --slide)")
    ap.add_argument("--target", action="store_true",
                    help="also draw the scoring curve, to <stem>_with_target.png")
    ap.add_argument("--slide", action="store_true",
                    help=f"{SLIDE_W_IN:g} in canvas for a slide, "
                         "to <stem>_slide.png")
    ap.add_argument("--loess", action="store_true",
                    help="draw the LOESS curve through the transects")
    ap.add_argument("--x", choices=("transect", "metres"), default="transect",
                    help="x axis: transect order (default) or along-coast metres")
    ap.add_argument("--no-domains", action="store_true",
                    help="transects only: no domain bands, numbers or means")
    a = ap.parse_args(argv)
    if a.gis is None:
        a.gis = list(SLIDE_GIS if a.slide else DEFAULT_GIS)
    apply_style()
    for p in figure(a.window, a.gis[0], a.gis[1], with_target=a.target,
                   slide=a.slide, with_loess=a.loess, x=a.x,
                   domains_shown=not a.no_domains):
        print(f"wrote    {p.relative_to(_REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
