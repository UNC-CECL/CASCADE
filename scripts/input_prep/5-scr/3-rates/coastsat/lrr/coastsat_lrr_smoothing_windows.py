"""
coastsat_lrr_smoothing_windows.py
==============================================================================
The three LOESS windows overlaid on ONE LRR rate field, in the house style
(Hannah, 2026-09-21). The rate is the thing smoothed here -- this is the field
itself, not a projection into metres and not a model comparison.

    coastsat/lrr/<w>/smoothing_windows_<w>.png

WHAT IT SHOWS
    The unsmoothed per-domain mean rate as the palest line, then the LOESS
    curve at 3, 5 and 10 domains (1.5, 2.5, 5.0 km) on a light-to-dark ramp,
    darkest being the window every run is actually graded at. Village bands,
    groin and piers, the hatched shoal boxes and the model-input beach fills
    come from the coastsat_lrr_windows panel, so this figure reads directly
    against lrr_<w>.png beside it -- same y bound, same marks.

WHY THE SIGN FILL IS NOT HERE
    lrr_<w>.png colours the rate blue seaward / red landward. Four curves
    share this panel, so that pair is not available: the smoothing width is an
    ORDERED variable and takes a sequential ramp instead, the one
    smoothing_scale.py uses, anchored on the house shoreline blue. Sign is
    read off the zero line, which is drawn.

THE SPLICE
    Every curve keeps the raw domain means over GIS 1-10 (the Oregon Inlet
    boundary treatment, cascade_pipeline.coastsat_loess.DEFAULT_LOESS), so all
    four are identical there by construction -- the same splice the scoring
    target is built through.

USAGE
    python scripts/input_prep/5-scr/3-rates/coastsat/lrr/coastsat_lrr_smoothing_windows.py
    python scripts/input_prep/5-scr/3-rates/coastsat/lrr/coastsat_lrr_smoothing_windows.py --window 1996_2024
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
from cascade_pipeline.coastsat_loess import (  # noqa: E402
    DEFAULT_LOESS, spliced_loess_series,
)
from site_layer.hat_figure_style import (  # noqa: E402
    DOMAIN_AXIS_LABEL, INK_MUTED, SMOOTH_RAMP, apply_style, caption, figsize,
    save,
)
from site_layer.hat_observed_rates import COASTSAT_LRR_ROOT, windows  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM  # noqa: E402

N = cw.N_DOMAINS
DEFAULT_WINDOW = "1996_2024"
# Domain units. 0 is the unsmoothed domain means; 10 is the grading window
# (cascade_pipeline.hindcast TARGET_WINDOW). 3 is roughly the alongshore
# decorrelation scale of the domain-mean rate, 1.5 km -- the set
# smoothing_scale/ sweeps (Hannah, 2026-09-21).
SMOOTH_WINDOWS = (0, 3, 5, 10)
GRADED_WINDOW = 10
# The two shoal-fronted peaks of the full-period field, the ones the widest
# window flattens. Read off the figure in smoothing_scale/PROVENANCE.md; the
# caption states the share of removed signal that actually falls in them
# rather than asserting the attribution.
PEAK_SPANS = ((29, 35), (65, 72))
LW_RAW = 0.8
LW_SMOOTH = 1.5
# The transect cloud the curves are fitted to (Hannah, 2026-09-22): the same
# size and alpha as the dots on lrr_<w>.png, but ONE neutral grey rather than
# the sign pair -- sign is the other figure's variable, smoothing width is
# this one's, and a second colour scale on the same panel reads as a third.
DOT_C = INK_MUTED
DOT_ALPHA = 0.35
DOT_S = 3.0


def km_of(window_domains):
    """A width in domain units as kilometres."""
    return window_domains * DOM.domain_spacing_m / 1000.0


def win_label(w):
    """One smoothing width, as the legend says it: the physical width alone.

    Which width the runs are graded at is a caption matter, not a legend one
    -- spelling it here ran the four entries past the figure edge. The domain
    count went the same way when the transect cloud took a fifth entry
    (2026-09-22); the caption gives both units for every width.
    """
    if not w:
        return "Unsmoothed domain means"
    return f"LOESS {km_of(w):g} km"


def transects(stem):
    """(frame, domain ids, along-coast metres, rate) for one LRR window.

    along-coast metres follow the convention every target build uses: each
    domain's transects spread evenly across its 500 m band, ordered within the
    domain by transect_id.
    """
    t = pd.read_csv(COASTSAT_LRR_ROOT / stem / "transect_lrr_full.csv")
    t = t[t["domain_number"].between(DOM.first_gis_id, DOM.last_gis_id)].copy()
    t["domain_number"] = t["domain_number"].astype(int)
    t = t.sort_values(["domain_number", "transect_id"]).reset_index(drop=True)
    rank = t.groupby("domain_number").cumcount()
    n = t.groupby("domain_number")["domain_number"].transform("count")
    sp = DOM.domain_spacing_m
    along = ((t["domain_number"] - DOM.first_gis_id) * sp
             + (rank + 0.5) * (sp / n)).to_numpy(float)
    return t, t["domain_number"].to_numpy(int), along, t["lrr_m_yr"].to_numpy(float)


def domain_x(along):
    """Along-coast metres on the GIS-domain axis the curves are drawn on."""
    return along / DOM.domain_spacing_m + DOM.first_gis_id - 0.5


def dots(ax, x, y, half):
    """The individual transect rates under the curves. Returns how many fell
    outside the shared y bound, drawn as open markers at the edge as they are
    on lrr_<w>.png."""
    ok = np.isfinite(y)
    inside = ok & (np.abs(y) <= half)
    ax.scatter(x[inside], y[inside], s=DOT_S, c=DOT_C, alpha=DOT_ALPHA,
               linewidths=0, zorder=3.5)
    out = ok & ~inside
    if out.any():
        ax.scatter(x[out], np.clip(y[out], -half * 0.985, half * 0.985), s=9,
                   facecolors="none", edgecolors=DOT_C, linewidths=0.6,
                   zorder=3.5)
    return int(out.sum())


def series(stem, smooth_windows=SMOOTH_WINDOWS):
    """{width: per-domain Series} plus {width: lowess frac}, spliced."""
    _, ids, along, rate = transects(stem)
    out, fracs = {}, {}
    for w in smooth_windows:
        out[w], fracs[w] = spliced_loess_series(ids, along, rate, w)
    return out, fracs


def structure(stem, vals, smooth_windows=SMOOTH_WINDOWS):
    """What each width takes out of THIS field, for the caption.

    Everything is in m/yr of alongshore structure removed, never a share of
    variance: the domain-mean variance is dominated by the long-wavelength
    swings, so a wiggle that is plainly visible on the figure reads as a few
    per cent of it and the percentage badly undersells the effect (Hannah,
    2026-09-21). These reproduce the target row of
    output/comparisons/target_comparison/smoothing_scale/tables/
    target_structure.csv exactly, and are recomputed here so the caption
    cannot drift from the table.
    """
    t, _, _, _ = transects(stem)
    dmean = t.groupby("domain_number")["lrr_m_yr"].mean()
    resid = {w: dmean.reindex(vals[w].index) - vals[w]
             for w in smooth_windows if w}
    widest = resid[smooth_windows[-1]]
    # Does the widest window take its structure out of the shoal-fronted
    # peaks, or evenly along the island? Squared removed signal per domain,
    # peaks against the rest, over the unspliced reach only.
    free = widest.loc[DEFAULT_LOESS.skip_southern_domains + 1:]
    peaks = [g for lo, hi in PEAK_SPANS for g in range(lo, hi + 1)]
    ss = free ** 2
    inpk = float(ss[ss.index.isin(peaks)].sum())
    return dict(
        sd_domain_mean=float(dmean.std(ddof=1)),
        sd_within=float((t["lrr_m_yr"] - t["domain_number"].map(dmean)).std(ddof=1)),
        median_unc=float(t["unc_m_yr"].median()),
        sd_removed={w: float(r.std(ddof=1)) for w, r in resid.items()},
        peak_share=100.0 * inpk / float(ss.sum()),
        n_peak=len(peaks), n_free=int(len(ss)),
    )


def lrr_half():
    """The y bound of every LRR window figure, so this one matches the figure
    already beside it (rates_figures.lrr_figures)."""
    frames = [pd.read_csv(COASTSAT_LRR_ROOT / f"{s}_{e}" / "domain_lrr_summary.csv")
              for s, e in windows()]
    return cw.shared_bounds(
        [rf._frame(f.assign(domain_number=f.domain_number.astype(int)), "mean_lrr")
         for f in frames])


def figure(stem, half=None, smooth_windows=SMOOTH_WINDOWS):
    """Draw one window's smoothing family. Returns the paths written."""
    start, end = (int(v) for v in stem.split("_"))
    half = lrr_half() if half is None else half
    vals, fracs = series(stem, smooth_windows)
    skip = DEFAULT_LOESS.skip_southern_domains

    fig, ax = plt.subplots(figsize=figsize("double", aspect=0.40),
                           constrained_layout=True)
    # mean_lrr all-NaN draws the frame, grid, village bands and structures with
    # no sign fill: the ramp below carries the ordered variable instead.
    frame = pd.DataFrame({"domain_number": np.arange(1, N + 1),
                          "mean_lrr": np.nan, "std_lrr": 0.0})
    cw.draw_panel(ax, frame, half, std=False)
    cw.draw_shoals(ax, label=True)
    fills = cw.fills_in(start, end)
    if fills:
        cw.draw_fills(ax, fills, half)

    # The cloud first, so every curve reads on top of the data it is fitted to.
    _, _, along, rate = transects(stem)
    n_dots_out = dots(ax, domain_x(along), rate, half)

    n_out = 0
    for i, w in enumerate(smooth_windows):
        s = vals[w]
        x = s.index.to_numpy(float)
        y = s.to_numpy(float)
        n_out += int(np.sum(np.isfinite(y) & (np.abs(y) > half)))
        colour = SMOOTH_RAMP[i % len(SMOOTH_RAMP)]
        if w:
            ax.plot(x, y, color=colour, lw=LW_SMOOTH, zorder=10 + i,
                    solid_capstyle="round")
        else:
            ax.plot(x, y, color=colour, lw=LW_RAW, marker="o", ms=1.8,
                    zorder=10, solid_capstyle="round")

    ax.yaxis.set_major_locator(MultipleLocator(cw.Y_TICK_M))
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    ax.set_ylabel(cw.Y_LABEL)

    h = [Line2D([], [], color=DOT_C, marker="o", ms=2.2, lw=0, alpha=0.7),
         Line2D([], [], color=SMOOTH_RAMP[0], lw=LW_RAW, marker="o", ms=2.2)]
    h += [Line2D([], [], color=SMOOTH_RAMP[i % len(SMOOTH_RAMP)], lw=LW_SMOOTH)
          for i, w in enumerate(smooth_windows) if w]
    labels = ["Individual transects"] + [win_label(w) for w in smooth_windows]
    fig.legend(h, labels, loc="outside lower center", ncol=len(h), frameon=False,
               handler_map={tuple: HandlerTuple(ndivide=None, pad=0.3)})

    widths = ", ".join(f"{km_of(w):g} km ({w} domains)"
                       for w in smooth_windows[1:-1])
    fr = "; ".join(f"{km_of(w):g} km frac {fracs[w]:.3f}"
                   for w in smooth_windows if w)
    st = structure(stem, vals, smooth_windows)
    removed = "; ".join(f"{km_of(w):g} km {st['sd_removed'][w]:.3f}"
                        for w in smooth_windows if w)
    peak_txt = ", ".join(f"GIS {lo}–{hi}" for lo, hi in PEAK_SPANS)
    caption(fig, (
        f"Observed shoreline change rate by GIS domain (1 at Cape Point, 90 at Pea "
        f"Island), {start}–{end}, at the three alongshore smoothing widths. For each "
        "CoastSat transect the rate is the ordinary-least-squares slope of shoreline "
        f"position against date over the calendar window (1 January {start} to 31 "
        f"December {end}), seaward positive. The grey dots are those individual "
        "transect rates, the cloud every curve here is fitted to"
        + (f" ({n_dots_out} beyond the axis, drawn as open circles at its edge)"
           if n_dots_out else "")
        + ". The palest line with markers is the mean "
        "of the ~10 transects in each 500 m domain, unsmoothed; the three heavier "
        f"curves are a LOESS fitted to the transects at {widths} and "
        f"{km_of(smooth_windows[-1]):g} km ({smooth_windows[-1]} domains), light "
        "to dark, the darkest being the "
        f"{smooth_windows[-1]}-domain window every run is graded at "
        f"({fr}). Every curve keeps the raw domain means over GIS 1–{skip} — the "
        "boundary treatment at Oregon Inlet — so all four are identical there by "
        "construction. Each window takes this much alongshore structure out of the "
        f"field, as the standard deviation of what it removes: {removed} m/yr, "
        f"against a domain-mean standard deviation of {st['sd_domain_mean']:.3f} m/yr. "
        "Widening the window is not denoising this field: the scatter between "
        f"transects inside one domain is only {st['sd_within']:.3f} m/yr and is larger "
        f"than the median per-transect regression uncertainty ({st['median_unc']:.3f} "
        "m/yr), so it is not clean estimation noise, and the domain-mean rate "
        f"decorrelates alongshore at about 1.5 km, three times narrower than the "
        f"{km_of(smooth_windows[-1]):g} km window. What that window removes is "
        f"concentrated in the two shoal-fronted peaks ({peak_txt}), which hold "
        f"{st['peak_share']:.0f}% of the removed signal by squared amplitude on "
        f"{st['n_peak']} of the {st['n_free']} smoothed domains. "
        + rf._marks_clause(start, end)
        + f" The y axis is ±{half:g} m/yr, the bound of every window figure."
        + (" This window is context: no model run is graded against it."
           if (start, end) == (1996, 2024) else "")))
    out = save(fig, COASTSAT_LRR_ROOT / stem / f"smoothing_windows_{stem}")
    plt.close(fig)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--window", default=DEFAULT_WINDOW,
                    help="LRR window stem, e.g. 1996_2024, or 'all'")
    a = ap.parse_args(argv)
    apply_style()
    half = lrr_half()
    stems = ([f"{s}_{e}" for s, e in windows()] if a.window == "all"
             else [a.window])
    for stem in stems:
        for p in figure(stem, half=half):
            print(f"wrote    {p.relative_to(_REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
