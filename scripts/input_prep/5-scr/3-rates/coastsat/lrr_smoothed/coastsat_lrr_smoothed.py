"""
The LRR field at two LOWESS widths, both on one figure per window, with the table behind it.

    python scripts/input_prep/5-scr/3-rates/coastsat/lrr_smoothed/coastsat_lrr_smoothed.py
    python scripts/input_prep/5-scr/3-rates/coastsat/lrr_smoothed/coastsat_lrr_smoothed.py --windows 1996_2015 --widths 5 7

Reads the window's lrr tables (coastsat_domain_lrr.py) and writes to 3-rates/coastsat/lrr_smoothed/<window>/.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-02
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
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import MultipleLocator  # noqa: E402

import coastsat_lrr_smoothing_windows as sw  # noqa: E402
import coastsat_lrr_windows as cw  # noqa: E402
import rates_figures as rf  # noqa: E402
from cascade_pipeline.coastsat_lowess import (  # noqa: E402
    DEFAULT_LOWESS, spliced_lowess_series,
)
from site_layer.hat_figure_style import (  # noqa: E402
    C, DOMAIN_AXIS_LABEL, SMOOTH_RAMP, _title, apply_style, caption, figsize, save,
    support_dir,
)
from site_layer.hat_observed_rates import COASTSAT_LRR_SMOOTHED_ROOT  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
# The candidate windows and their reference (Hannah, 2026-10-02)
DEFAULT_WINDOWS = ("1996_2015", "2010_2026", "1996_2026")
# Domain units; 7 is the group's range and the graded width, 5 the narrower test
DEFAULT_WIDTHS = (5, 7)
# Raw means palest, then the widths light to dark
COLOURS = {0: SMOOTH_RAMP[0], "narrow": SMOOTH_RAMP[2], "wide": SMOOTH_RAMP[3]}
# The long-term rate drawn on each candidate, at the wider width (Hannah, 2026-10-02)
DEFAULT_REFERENCE = "1996_2026"
REF_C = C["BASE"]
REF_LW = 1.3
REF_DASH = (0, (4, 2.2))
# The 95% bands, faint enough that the curves stay on top
BAND_ALPHA = 0.20
# The 2021 step test: the same start, the record cut before the step (Hannah, 2026-10-02)
STEP_WINDOWS = ("2010_2020", "2010_2026")
STEP_COLOURS = {"2010_2020": C["ACCENT"], "2010_2026": SMOOTH_RAMP[3]}
# The two candidates side by side, in the colours of the lrr overlay figure (Hannah, 2026-10-02)
CANDIDATE_PAIR = ("1996_2015", "2010_2026")
CANDIDATE_COLOURS = {"1996_2015": "0.15", "2010_2026": C["ACCENT"]}
COLOUR_NAMES = {C["ACCENT"]: "purple", SMOOTH_RAMP[3]: "dark blue", "0.15": "black"}
# The narrower width as a tint of the window's colour, in the widths overlay
NARROW_TINT = 0.55
# -----------------------------------------------------------------------------


# The colour of each curve, raw first
def colour_of(i, widths):
    return COLOURS["narrow"] if widths[i] == min(widths) else COLOURS["wide"]


# Raw domain means and each width's LOWESS as one table, with the lowess fracs
def smoothed_table(stem, widths):
    _, ids, along, rate = sw.transects(stem)
    raw, _ = spliced_lowess_series(ids, along, rate, 0)
    table = pd.DataFrame({"mean_lrr": raw.round(3)})
    fracs = {}
    for w in widths:
        s, fracs[w] = spliced_lowess_series(ids, along, rate, w)
        table[f"lowess{w}"] = s.round(3)
    a, b = widths
    table[f"lowess{b}_minus_lowess{a}"] = (table[f"lowess{b}"]
                                           - table[f"lowess{a}"]).round(3)
    return table.reset_index(), fracs


# What separates the two widths, for the caption
def difference_stats(table, widths):
    a, b = widths
    d = table.set_index("domain_number")[f"lowess{b}_minus_lowess{a}"]
    d = d.loc[DEFAULT_LOWESS.skip_southern_domains + 1:].dropna()
    return dict(sd=float(d.std(ddof=1)), max_abs=float(d.abs().max()),
                at=int(d.abs().idxmax()), n_over_half=int((d.abs() > 0.5).sum()),
                n=int(len(d)))


# A window's rate and its 95% half-width, both smoothed at one width
def curve_and_band(stem, width):
    t, ids, along, rate = sw.transects(stem)
    s, _ = spliced_lowess_series(ids, along, rate, width)
    u, _ = spliced_lowess_series(ids, along, t["unc_m_yr"].to_numpy(float),
                                 width)
    return s.round(3), u.round(3)


# How far one curve sits from another north of the splice, and where their bands do not overlap
def offset_stats(s, u, s_ref, u_ref):
    lo = DEFAULT_LOWESS.skip_southern_domains + 1
    d = (s - s_ref).loc[lo:]
    gap = d.abs() > (u + u_ref).loc[lo:]
    ok = d.notna()
    return dict(mean=float(d[ok].mean()), rms=float(np.sqrt((d[ok] ** 2).mean())),
                n=int(ok.sum()), n_apart=int(gap[ok].sum()),
                apart=[int(g) for g in gap[ok & gap].index])


# Domain ids as compact runs, "GIS 7–19, 63–66"
def runs_text(ids):
    if not ids:
        return "none"
    runs, start, prev = [], ids[0], ids[0]
    for g in ids[1:] + [None]:
        if g is not None and g == prev + 1:
            prev = g
            continue
        runs.append(f"{start}" if start == prev else f"{start}–{prev}")
        if g is not None:
            start = prev = g
    return "GIS " + ", ".join(runs)


# Above or below, with the size
def side(v):
    return f"{abs(v):.2f} m/yr {'above' if v >= 0 else 'below'}"


# A curve and its band
def draw_curve(ax, x, s, u, colour, lw, zorder, ls="-"):
    ax.fill_between(x, s - u, s + u, color=colour, alpha=BAND_ALPHA, lw=0,
                    zorder=zorder - 0.5)
    ax.plot(x, s, color=colour, lw=lw, ls=ls, zorder=zorder,
            solid_capstyle="round")


# One window: the transect cloud, the raw domain means, both widths and the reference
def figure(stem, half, widths, ref=None):
    start, end = (int(v) for v in stem.split("_"))
    table, fracs = smoothed_table(stem, widths)
    st = difference_stats(table, widths)
    skip = DEFAULT_LOWESS.skip_southern_domains
    ref = None if ref == stem else ref
    a, b = min(widths), max(widths)
    _, unc_b = curve_and_band(stem, b)
    table[f"unc95_lowess{b}"] = unc_b.to_numpy()
    if ref:
        r0, r1 = (int(v) for v in ref.split("_"))
        ref_s, ref_u = curve_and_band(ref, b)
        table[f"ref_lowess{b}_{ref}"] = ref_s.to_numpy()
        table[f"ref_unc95_lowess{b}_{ref}"] = ref_u.to_numpy()
        rs = offset_stats(table.set_index("domain_number")[f"lowess{b}"], unc_b,
                          ref_s, ref_u)

    fig, ax = plt.subplots(figsize=figsize("double", aspect=0.40),
                           constrained_layout=True)
    frame = pd.DataFrame({"domain_number": np.arange(1, cw.N_DOMAINS + 1),
                          "mean_lrr": np.nan, "std_lrr": 0.0})
    cw.draw_panel(ax, frame, half, std=False)
    cw.draw_shoals(ax, label=True)
    fills = cw.fills_in(start, end)
    if fills:
        cw.draw_fills(ax, fills, half)

    _, _, along, rate = sw.transects(stem)
    n_dots_out = sw.dots(ax, sw.domain_x(along), rate, half)

    x = table["domain_number"].to_numpy(float)
    ax.plot(x, table["mean_lrr"], color=COLOURS[0], lw=sw.LW_RAW, marker="o",
            ms=1.8, zorder=10, solid_capstyle="round")
    for i, w in enumerate(widths):
        if w == b:
            draw_curve(ax, x, table[f"lowess{w}"], table[f"unc95_lowess{w}"],
                       colour_of(i, widths), sw.LW_SMOOTH, 11 + i)
        else:
            ax.plot(x, table[f"lowess{w}"], color=colour_of(i, widths),
                    lw=sw.LW_SMOOTH, zorder=11 + i, solid_capstyle="round")
    if ref:
        draw_curve(ax, x, table[f"ref_lowess{b}_{ref}"],
                   table[f"ref_unc95_lowess{b}_{ref}"], REF_C, REF_LW,
                   11 + len(widths), ls=REF_DASH)

    ax.yaxis.set_major_locator(MultipleLocator(cw.Y_TICK_M))
    ax.set_title(f"Shoreline change rate, {start}–{end}: LOWESS smoothing at "
                 f"{sw.km_of(widths[0]):g} km and {sw.km_of(widths[1]):g} km",
                 loc="center", pad=cw.TITLE_PAD_FILLS if fills else None)
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    ax.set_ylabel(cw.Y_LABEL)

    h = [Line2D([], [], color=sw.DOT_C, marker="o", ms=2.2, lw=0, alpha=0.7),
         Line2D([], [], color=COLOURS[0], lw=sw.LW_RAW, marker="o", ms=2.2)]
    h += [Line2D([], [], color=colour_of(i, widths), lw=sw.LW_SMOOTH)
          for i in range(len(widths))]
    labels = (["Individual transects", sw.win_label(0)]
              + [f"LOWESS {sw.km_of(w):g} km ({w} domains)" for w in widths])
    if ref:
        h.append(Line2D([], [], color=REF_C, lw=REF_LW, ls=REF_DASH))
        labels.append(f"{r0}–{r1}, LOWESS {sw.km_of(max(widths)):g} km")
    # Four fit on one row of a double-column figure; five go to two
    fig.legend(h, labels, loc="outside lower center",
               ncol=len(h) if len(h) <= 4 else 3, frameon=False)

    caption(fig, (
        f"Observed shoreline change rate by GIS domain (1 at Cape Point, 90 at "
        f"Pea Island), {start}–{end}, smoothed alongshore at two widths. For "
        "each CoastSat transect the rate is the ordinary-least-squares slope "
        "of shoreline position against date over the calendar window (1 "
        f"January {start} to 31 December {end}), seaward positive; the grey "
        "dots are those transect rates"
        + (f" ({n_dots_out} beyond the axis, drawn as open circles at its edge)"
           if n_dots_out else "")
        + ". The palest line with markers is the mean of the ~10 transects in "
        f"each 500 m domain, unsmoothed. The two heavier curves are a LOWESS "
        f"fitted to the transects at {sw.km_of(a):g} km ({a} domains, frac "
        f"{fracs[a]:.3f}) and {sw.km_of(b):g} km ({b} domains, frac "
        f"{fracs[b]:.3f}), the darker being the {b}-domain width every run is "
        f"graded at. Both keep the raw domain means over GIS 1–{skip}, the "
        "boundary treatment at Oregon Inlet, so they are identical there by "
        f"construction. North of GIS {skip} the {b}-domain curve differs from "
        f"the {a}-domain one by {st['sd']:.3f} m/yr (standard deviation over "
        f"{st['n']} domains), at most {st['max_abs']:.2f} m/yr at GIS "
        f"{st['at']}, and by more than 0.5 m/yr at {st['n_over_half']} "
        "domains. "
        + f"The shaded band on the {b}-domain curve is the 95% confidence "
        "interval of the transect fits, smoothed the same way (the half-widths "
        "are averaged, not combined as independent errors, because neighbouring "
        "transects are not independent). "
        + (f"The dashed grey line and band are the long-term {r0}–{r1} rate, "
           f"fitted and smoothed the same way. North of GIS {skip} this "
           f"window's {b}-domain curve sits {side(rs['mean'])} it on average "
           f"(root mean square difference {rs['rms']:.2f} m/yr), and the two "
           f"bands do not overlap at {rs['n_apart']} of {rs['n']} domains "
           f"({runs_text(rs['apart'])}). "
           if ref else "")
        + rf._marks_clause(start, end)
        + f" The y axis is ±{half:g} m/yr, shared by every figure of the "
        "candidate windows."))

    out_dir = COASTSAT_LRR_SMOOTHED_ROOT / stem
    written = save(fig, out_dir / f"lrr_lowess{a}_vs_lowess{b}_{stem}")
    plt.close(fig)
    csv = support_dir(out_dir) / f"lrr_lowess{a}_vs_lowess{b}_{stem}.csv"
    table.to_csv(csv, index=False)
    return written + [csv], st, (rs if ref else None)


# Two windows over the reference, one panel per width; `about` says why the pair is drawn
def pair_figure(half, widths, ref, windows, colours, about, out_dir):
    short, long_ = windows
    (s0, s1), (l0, l1) = ((int(v) for v in w.split("_")) for w in windows)
    r0, r1 = (int(v) for v in ref.split("_"))
    skip = DEFAULT_LOWESS.skip_southern_domains
    fills = cw.fills_in(min(s0, l0), max(s1, l1))
    curves, st, pair = {}, {}, {}
    for width in widths:
        curves[width] = {w: curve_and_band(w, width) for w in (*windows, ref)}
        c = curves[width]
        st[width] = {w: offset_stats(*c[w], *c[ref]) for w in windows}
        pair[width] = offset_stats(*c[long_], *c[short])

    n = len(widths)
    fig, axes = plt.subplots(
        n, 1, sharex=True, sharey=True, squeeze=False, constrained_layout=True,
        figsize=(figsize("double", aspect=0.40) if n == 1
                 else figsize("double", height=5.4)))
    axes = axes[:, 0]
    frame = pd.DataFrame({"domain_number": np.arange(1, cw.N_DOMAINS + 1),
                          "mean_lrr": np.nan, "std_lrr": 0.0})
    for k, (ax, width) in enumerate(zip(axes, widths)):
        cw.draw_panel(ax, frame, half, std=False, label=(k == 0))
        cw.draw_shoals(ax, label=(k == 0))
        if fills and k == 0:
            cw.draw_fills(ax, fills, half)
        c = curves[width]
        x = c[ref][0].index.to_numpy(float)
        draw_curve(ax, x, *c[ref], REF_C, REF_LW, 10, ls=REF_DASH)
        for i, w in enumerate(windows):
            draw_curve(ax, x, *c[w], colours[w], sw.LW_SMOOTH, 11 + i)
        ax.yaxis.set_major_locator(MultipleLocator(cw.Y_TICK_M))
        if n > 1:
            _title(ax, k, f"LOWESS {sw.km_of(width):g} km ({width} domains)")
    axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)
    head = (f"Shoreline change rate, {s0}–{s1} and {l0}–{l1}, against the "
            f"{r0}–{r1} rate")
    if n == 1:
        axes[0].set_title(head + f" (LOWESS {sw.km_of(widths[0]):g} km)",
                          loc="center", pad=cw.TITLE_PAD_FILLS if fills else None)
        axes[0].set_ylabel(cw.Y_LABEL)
    else:
        fig.suptitle(head, fontsize=10)
        fig.supylabel(cw.Y_LABEL, fontsize=9)
    h = [Line2D([], [], color=colours[w], lw=sw.LW_SMOOTH) for w in windows]
    h.append(Line2D([], [], color=REF_C, lw=REF_LW, ls=REF_DASH))
    fig.legend(h, [f"{s0}–{s1}", f"{l0}–{l1}", f"{r0}–{r1}"],
               loc="outside lower center", ncol=3, frameon=False)

    def stats_txt(width):
        a, b, p_ = st[width][short], st[width][long_], pair[width]
        return (f"{s0}–{s1} sits {side(a['mean'])} the long-term rate on average "
                f"(bands apart at {a['n_apart']} of {a['n']} domains) and "
                f"{l0}–{l1} sits {side(b['mean'])} it (bands apart at "
                f"{b['n_apart']}); {l0}–{l1} sits {side(p_['mean'])} {s0}–{s1} "
                f"(root mean square difference {p_['rms']:.2f} m/yr, bands apart "
                f"at {p_['n_apart']} domains: {runs_text(p_['apart'])})")

    widths_txt = " and ".join(f"{sw.km_of(w):g} km ({w} domains)" for w in widths)
    if n == 1:
        stats = f"North of GIS {skip}, " + stats_txt(widths[0]) + ". "
    else:
        stats = f"North of GIS {skip}: " + "; ".join(
            f"({chr(97 + k)}) at {w} domains, {stats_txt(w)}"
            for k, w in enumerate(widths)) + ". "
    caption(fig, (
        f"Observed shoreline change rate by GIS domain (1 at Cape Point, 90 at "
        f"Pea Island) for {s0}–{s1} ({COLOUR_NAMES[colours[short]]}) and "
        f"{l0}–{l1} ({COLOUR_NAMES[colours[long_]]}), {about}, over the "
        f"long-term {r0}–{r1} rate (dashed grey). Each curve is "
        "the per-transect CoastSat linear regression rate over the calendar "
        f"window, LOWESS-smoothed at {widths_txt}"
        + (", one panel each," if n > 1 else "")
        + f" with the raw domain means kept over GIS 1–{skip}; each band is the "
        "95% confidence interval of the transect fits, smoothed the same way. "
        + stats
        + (f"What separates the two windows comes from {s1 + 1}–{l1}, the years "
           "they do not share. " if s0 == l0 else
           f"The two share {max(s0, l0)}–{min(s1, l1)}. ")
        + f"Black bars above the frame mark the fills inside "
        f"{min(s0, l0)}–{max(s1, l1)} at the hindcast footprint. "
        + rf._marks_clause(min(s0, l0), max(s1, l1))
        + f" The y axis is ±{half:g} m/yr, shared by every figure of the "
        "candidate windows."))

    wtag = "_and_".join(f"lowess{w}" for w in widths)
    stem = f"lrr_{wtag}_{short}_and_{long_}_vs_{ref}"
    written = save(fig, out_dir / stem)
    plt.close(fig)
    table = pd.DataFrame({"domain_number": x.astype(int)})
    for width in widths:
        for w in (*windows, ref):
            table[f"lowess{width}_{w}"] = curves[width][w][0].to_numpy()
            table[f"unc95_lowess{width}_{w}"] = curves[width][w][1].to_numpy()
    csv = support_dir(out_dir) / f"{stem}.csv"
    table.to_csv(csv, index=False)
    return written + [csv], st, pair


# A hex or grey-level colour mixed toward white
def tint(colour, amount):
    r, g, b = matplotlib.colors.to_rgb(colour)
    return tuple(c + (1 - c) * amount for c in (r, g, b))


# Both candidate windows at both widths on one panel (Hannah, 2026-10-02)
def widths_overlay_figure(half, widths, windows, colours, out_dir):
    a, b = widths
    skip = DEFAULT_LOWESS.skip_southern_domains
    curves = {w: {k: curve_and_band(w, k)[0] for k in widths} for w in windows}
    fig, ax = plt.subplots(figsize=figsize("double", aspect=0.40),
                           constrained_layout=True)
    frame = pd.DataFrame({"domain_number": np.arange(1, cw.N_DOMAINS + 1),
                          "mean_lrr": np.nan, "std_lrr": 0.0})
    cw.draw_panel(ax, frame, half, std=False)
    cw.draw_shoals(ax, label=True)
    stats = {}
    h, labels = [], []
    for k, w in enumerate(windows):
        x = curves[w][a].index.to_numpy(float)
        ax.plot(x, curves[w][a], color=tint(colours[w], NARROW_TINT),
                lw=sw.LW_SMOOTH, zorder=11 + 2 * k, solid_capstyle="round")
        ax.plot(x, curves[w][b], color=colours[w], lw=sw.LW_SMOOTH * 0.8,
                zorder=12 + 2 * k, solid_capstyle="round")
        d = (curves[w][b] - curves[w][a]).loc[skip + 1:]
        stats[w] = dict(sd=float(d.std(ddof=1)), max_abs=float(d.abs().max()),
                        at=int(d.abs().idxmax()),
                        n_half=int((d.abs() > 0.5).sum()), n=int(d.notna().sum()))
        h += [Line2D([], [], color=tint(colours[w], NARROW_TINT), lw=sw.LW_SMOOTH),
              Line2D([], [], color=colours[w], lw=sw.LW_SMOOTH * 0.8)]
        labels += [f"{w.replace('_', '–')}, {sw.km_of(a):g} km ({a} domains)",
                   f"{w.replace('_', '–')}, {sw.km_of(b):g} km ({b} domains)"]
    ax.yaxis.set_major_locator(MultipleLocator(cw.Y_TICK_M))
    ax.set_title(f"Shoreline change rate, " + " and ".join(
        w.replace("_", "–") for w in windows)
        + f", at LOWESS {sw.km_of(a):g} km and {sw.km_of(b):g} km", loc="center")
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    ax.set_ylabel(cw.Y_LABEL)
    fig.legend(h, labels, loc="outside lower center", ncol=2, frameon=False)

    st_txt = "; ".join(
        f"{w.replace('_', '–')} standard deviation {v['sd']:.2f} m/yr, largest "
        f"{v['max_abs']:.2f} at GIS {v['at']}, more than 0.5 m/yr at "
        f"{v['n_half']} of {v['n']} domains" for w, v in stats.items())
    caption(fig, (
        "Observed shoreline change rate by GIS domain (1 at Cape Point, 90 at Pea "
        "Island) for the candidate windows "
        + " (black) and ".join(w.replace("_", "–") for w in windows)
        + " (purple), each LOWESS-smoothed alongshore at "
        f"{sw.km_of(a):g} km ({a} domains, the lighter line) and {sw.km_of(b):g} km "
        f"({b} domains, the darker line), all on one panel; the bands and the "
        "long-term rate of the side-by-side figure are left off so the widths "
        "can be told apart. The rates are per-transect CoastSat linear regression "
        "slopes over the calendar window; both widths keep the raw domain means "
        f"over GIS 1–{skip}, where they are identical by construction. The "
        f"{b}-domain curve minus the {a}-domain one, north of GIS {skip}: "
        f"{st_txt}. The y axis is ±{half:g} m/yr, shared by every figure of the "
        "candidate windows."))
    stem = (f"lrr_lowess{a}_over_lowess{b}_"
            + "_and_".join(windows))
    written = save(fig, out_dir / stem)
    plt.close(fig)
    table = pd.DataFrame({"domain_number": curves[windows[0]][a].index.astype(int)})
    for w in windows:
        for k in widths:
            table[f"lowess{k}_{w}"] = curves[w][k].to_numpy()
        table[f"lowess{b}_minus_lowess{a}_{w}"] = (curves[w][b]
                                                    - curves[w][a]).round(3).to_numpy()
    csv = support_dir(out_dir) / f"{stem}.csv"
    table.to_csv(csv, index=False)
    return written + [csv], stats


# Run: one figure and table per window, every window on one y axis
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 2)[1])
    ap.add_argument("--windows", nargs="+", default=list(DEFAULT_WINDOWS),
                    metavar="START_END")
    ap.add_argument("--widths", nargs=2, type=int, default=list(DEFAULT_WIDTHS),
                    metavar="DOMAINS", help="two LOWESS widths, narrower first")
    ap.add_argument("--reference", default=DEFAULT_REFERENCE, metavar="START_END",
                    help="long-term window drawn dashed at the wider width "
                         "('none' to leave it off)")
    a = ap.parse_args(argv)
    widths = tuple(sorted(a.widths))
    apply_style()
    half, _ = cw.candidate_half()
    ref = None if a.reference == "none" else a.reference
    for stem in a.windows:
        written, st, rs = figure(stem, half, widths, ref)
        print(f"{stem}  sd {st['sd']:.3f}  max {st['max_abs']:.2f} at GIS "
              f"{st['at']}  >0.5 at {st['n_over_half']}/{st['n']}"
              + (f"  vs {ref} mean {rs['mean']:+.2f} rms {rs['rms']:.2f} "
                 f"apart {rs['n_apart']}/{rs['n']} ({runs_text(rs['apart'])})"
                 if rs else ""))
        for p in written:
            print(f"wrote    {p.relative_to(_REPO)}")
    # The 2021-step check at the graded width; the candidate pair at both widths (Hannah, 2026-10-02)
    pairs = (
        (STEP_WINDOWS, STEP_COLOURS, (max(widths),), "the same start with the "
         "record cut before and after the 2021 shoreline step",
         COASTSAT_LRR_SMOOTHED_ROOT / "2010_2026"),
        (CANDIDATE_PAIR, CANDIDATE_COLOURS, widths, "the two candidate windows",
         COASTSAT_LRR_SMOOTHED_ROOT / ref if ref else None),
    )
    for windows, colours, pw, about, out_dir in (pairs if ref else ()):
        written, st, pair = pair_figure(half, pw, ref, windows, colours, about,
                                        out_dir)
        for width in pw:
            for w, v in st[width].items():
                print(f"[{width}] {w} vs {ref}  mean {v['mean']:+.2f}  "
                      f"apart {v['n_apart']}/{v['n']}")
            p_ = pair[width]
            print(f"[{width}] {windows[1]} vs {windows[0]}  mean {p_['mean']:+.2f}  "
                  f"rms {p_['rms']:.2f}  apart {p_['n_apart']} "
                  f"({runs_text(p_['apart'])})")
        for p in written:
            print(f"wrote    {p.relative_to(_REPO)}")
    if ref:
        written, stats = widths_overlay_figure(
            half, widths, CANDIDATE_PAIR, CANDIDATE_COLOURS,
            COASTSAT_LRR_SMOOTHED_ROOT / ref)
        for w, v in stats.items():
            print(f"[{widths[1]}-{widths[0]}] {w}  sd {v['sd']:.2f}  max "
                  f"{v['max_abs']:.2f} at GIS {v['at']}  >0.5 at {v['n_half']}/{v['n']}")
        for p in written:
            print(f"wrote    {p.relative_to(_REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
