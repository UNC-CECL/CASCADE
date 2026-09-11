r"""
HAT_method_comparison_figures.py
===============================================================================
The two setback methods against each other, and both against where NC-12
actually is.

  HAT_method_comparison_on_domains.png   both methods, each period on its own
                                         interiors
  HAT_method_vs_actual_road.png          both methods against the rasterized road

Written to road_offset/ itself rather than into either method's folder, because
neither figure belongs to one method.

ENCODING CHANGES HERE -- READ THIS FIRST
----------------------------------------
In the per-method figures, hue = YEAR. In these two, hue = METHOD and the year
is the panel:

    grey    the superseded method: the minimum road elevation minus the
            minimum dune elevation, taken independently per domain, against
            the same-year digitised dune line
    purple  the current method: per profile, measured landward from the dune
            start, then the domain median

That is the house BASE/ACCENT pair -- the input as it was against the change
under test -- and it leaves the vintage red/blue free. Because hue is NOT the
vintage in these two figures, the red pole is available, and it carries the
one thing that is a failure rather than a category: a roadway that drowns at
initialisation.

WHAT "ACTUAL ROAD" MEANS, AND WHAT IT DOES NOT
----------------------------------------------
The reference band is the RASTERIZED NC-12 mask, per alongshore profile, in the
same frame both methods are drawn in: `road_seaward_cell - interior_row0_cell`,
read from dunestart_offset/<year>/RoadOffset_<year>_profiles.csv. It is the road
as the model grid sees it, with all 50 profiles kept instead of collapsed.

It is NOT an independent check on the dune-start method. That method's setback
IS the median of this quantity, so the two agree by construction, differing only
by `int()` truncation and the negative floor. Read the band for two things it
does show:

  * how far the OLD method sits from the road actually burnt on the grid, which
    IS an independent comparison, because that method never saw this grid;
  * how much the road wanders WITHIN a domain -- the p10-p90 spread that any
    single scalar setback has to throw away, whichever method produced it.

REQUIREMENTS
------------
  numpy, pandas, matplotlib
===============================================================================
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.patheffects import withStroke

# =============================================================================
# SHARED MACHINERY
# =============================================================================
# The placement script owns load_years / load_interiors / place_road and the
# transcribed drown test. Importing it keeps ONE implementation: a second copy
# would drift, and a drifted method comparison looks like a result.

_HERE = Path(__file__).resolve().parent
_PLACEMENT = _HERE.parent / "1-produce" / "HAT_road_placement_on_domains.py"
if not _PLACEMENT.is_file():
    raise SystemExit(f"cannot find the placement module it shares code with:\n"
                     f"    {_PLACEMENT}\n"
                     f"If it moved, update _PLACEMENT here -- do NOT copy its "
                     f"place_road/drown test into this file.")
_SPEC = importlib.util.spec_from_file_location("hat_placement", _PLACEMENT)
P = importlib.util.module_from_spec(_SPEC)
sys.modules["hat_placement"] = P
_SPEC.loader.exec_module(P)

sys.path.insert(0, str(_HERE.parents[3]))
from hat_figure_style import (  # noqa: E402
    C, C_1984, DOMAIN_AXIS_LABEL, apply_style, caption, figsize, open_frame,
    save, spines_for_image, town_bands, _title)

apply_style()

ROADS_ROOT = P.ROADS_ROOT

# Cross-method output goes in its own folder, NOT at road_offset/ level and NOT
# inside either method's folder. The rule this satisfies is unchanged -- a
# legacy-vs-dune-start result belongs to neither method -- but the top level is
# for the forcing, its source and its inputs, and four loose comparison files
# sitting beside them read as though they were part of the product.
OUT_ROOT = ROADS_ROOT / "method_comparison"
YEARS = P.YEARS
DOMAINS = P.DOMAINS
CELL_SIZE_M = P.CELL_SIZE_M
DROWN_PCT = P.DROWN_PCT
SURFACE, WATER, NODATA = P.SURFACE, P.WATER, P.NODATA
INK_MUTED, INK_SECOND = P.INK_MUTED, P.INK_SECOND
INK = P.INK_SECOND

# hue = METHOD here, not year. See the header. BASE is the input as it stood,
# ACCENT the change under test. The rasterized road is the OBSERVATION both
# methods are measured against, which is what C["REF"] means -- it was the
# road ink for one draft, and a 42%-alpha near-black band is the same mid grey
# as BASE, so the reference and the superseded method read as one thing. A
# drowning roadway is the one failure state, and it takes the red pole, free
# in this figure precisely because hue is not the vintage here.
C_OLD, C_NEW = C["BASE"], C["ACCENT"]
C_ACTUAL = C["REF"]
C_DROWN = C_1984

METHOD_ORDER = ["old", "dunestart"]
METHOD_COLOUR = {"old": C_OLD, "dunestart": C_NEW}
# What the methods ARE, not what they were called while they were being built.
METHOD_LABEL = {
    "old": "setback from independent minima (superseded)",
    "dunestart": "setback from the dune start",
}
METHOD_TICK = {"old": "independent minima", "dunestart": "dune start"}

PROFILES_FMT = ("dunestart_offset/{year}/RoadOffset_{year}_profiles.csv")


# =============================================================================
# LOAD
# =============================================================================

def load_actual(year: int) -> dict:
    """
    The rasterized road per domain, kept as a spread rather than a scalar.

    seaward_p10/p50/p90 are percentiles ACROSS the domain's profiles of the
    road's seaward edge relative to interior row 0. width is the median
    measured road width, so the band drawn is a real footprint.
    """
    p = ROADS_ROOT / PROFILES_FMT.format(year=year)
    if not p.is_file():
        return {}
    df = pd.read_csv(p)
    out = {}
    for d, g in df.groupby("domain"):
        sb = g["setback_m"].to_numpy(dtype=float)
        out[int(d)] = dict(
            p10=float(np.percentile(sb, 10)),
            p50=float(np.median(sb)),
            p90=float(np.percentile(sb, 90)),
            width=float(np.median(g["road_width_m"])),
            n=int(len(sb)),
            spread=float(np.percentile(sb, 90) - np.percentile(sb, 10)),
        )
    return out


def load_placements(per: dict) -> dict:
    """{method: {year: placed}} using the placement script's own logic.

    `per` is the placement script's per-VINTAGE bundle, not one interiors dict.
    It used to be the latter, which placed the 1984 road on 2004-start
    interiors -- see the note at the top of HAT_road_placement_on_domains.py.
    """
    out = {}
    for name in METHOD_ORDER:
        spec = P.METHODS[name]
        per_year = {}
        for year in YEARS:
            sb = P.read_two_row(spec["root"] / spec["setback"].format(year=year))
            if sb:
                per_year[year] = P.place_road(per[year]["interiors"], sb)
        if per_year:
            out[name] = per_year
    return out


# =============================================================================
# SHARED DRAWING
# =============================================================================

def base_panel(ax, fig, shown, crop_rows, title, panel_index):
    ax.set_facecolor(NODATA)
    im = ax.imshow(np.ma.masked_invalid(shown), aspect="auto", origin="lower",
                   extent=[0.5, len(DOMAINS) + 0.5, -CELL_SIZE_M / 2,
                           crop_rows * CELL_SIZE_M - CELL_SIZE_M / 2],
                   cmap=P.LAND_CLASS_CMAP, norm=P.LAND_CLASS_NORM,
                   interpolation="nearest")
    cax = ax.inset_axes([1.012, 0.0, 0.014, 1.0])
    cb = fig.colorbar(im, cax=cax, spacing="uniform",
                      ticks=P.LAND_CLASS_BOUNDS[1:-1])
    cb.set_label("elevation (m MHW)")
    cb.outline.set_edgecolor(INK_MUTED)
    cb.outline.set_linewidth(0.6)

    ax.set_xlim(0.5, len(DOMAINS) + 0.5)
    if panel_index == 0:
        town_bands(ax, strip=0.075, shade=SURFACE)

    spines_for_image(ax)
    _title(ax, panel_index, title)
    ax.set_ylabel("m landward of\ninterior row 0")
    plt.setp(ax.get_xticklabels(), visible=False)


def step_xy(values: dict, key="start_m"):
    """One flat segment per domain, so a road is a road and not a polyline."""
    xs, ys = [], []
    for d in sorted(values):
        v = values[d] if not isinstance(values[d], dict) else values[d][key]
        xs += [d - 0.5, d + 0.5, np.nan]
        ys += [v, v, np.nan]
    return xs, ys


def draw_method_line(ax, placed, colour, lw=2.4, zorder=6, halo=True):
    xs, ys = step_xy(placed)
    kw = dict(path_effects=[withStroke(linewidth=lw + 1.8, foreground=SURFACE)]
              ) if halo else {}
    ax.plot(xs, ys, color=colour, lw=lw, solid_capstyle="butt", zorder=zorder,
            **kw)


# =============================================================================
# FIGURE 1 -- both methods, each period on its own interiors
# =============================================================================

def figure_methods(per, crop_rows, placements, out_png: Path):
    fig = plt.figure(figsize=figsize("double", height=8.0))
    gs = fig.add_gridspec(4, 1, height_ratios=[1.25, 1.25, 1.0, 0.62],
                          hspace=0.30, left=0.135, right=0.870,
                          top=0.958, bottom=0.090)
    axes = [fig.add_subplot(gs[0])]
    axes += [fig.add_subplot(gs[i], sharex=axes[0]) for i in (1, 2, 3)]
    ax84, ax04, ax_d, ax_s = axes

    for i, (ax, year) in enumerate(((ax84, YEARS[0]), (ax04, YEARS[1]))):
        base_panel(ax, fig, per[year]["shown"], crop_rows,
                   f"both methods in {year}", i)
        for name in METHOD_ORDER:
            if year in placements.get(name, {}):
                draw_method_line(ax, placements[name][year],
                                 METHOD_COLOUR[name],
                                 lw=2.8 if name == "old" else 2.0,
                                 zorder=6 if name == "old" else 7)

    # --- (c) how far apart, per domain -------------------------------------
    ax_d.axhline(0, color=INK_MUTED, lw=0.8, zorder=3)
    medians = {}
    for year, style in zip(YEARS, [(0, ()), (0, (5, 1.6))]):
        if not all(year in placements.get(m, {}) for m in METHOD_ORDER):
            continue
        o, n = (placements["old"][year], placements["dunestart"][year])
        common = sorted(set(o) & set(n))
        diff = [n[d]["setback_m"] - o[d]["setback_m"] for d in common]
        medians[year] = float(np.median(diff))
        ax_d.plot(common, diff, color=INK, lw=1.3, ls=style, zorder=5,
                  label=f"{year}")
    ax_d.set_ylabel("dune start − independent\nminima, setback (m)")
    ax_d.grid(axis="y")
    ax_d.set_axisbelow(True)
    town_bands(ax_d, label=False)
    open_frame(ax_d)
    ax_d.legend(loc="lower left", ncol=2, fontsize=7)
    _title(ax_d, 2, "difference between the two methods")
    plt.setp(ax_d.get_xticklabels(), visible=False)

    # --- (D) drown status, one row per method-year --------------------------
    rows = [(m, y) for m in METHOD_ORDER for y in YEARS
            if y in placements.get(m, {})]
    for k, (m, y) in enumerate(rows):
        pl = placements[m][y]
        bad = [d for d, p in pl.items() if p["drowned"]]
        ax_s.scatter(sorted(pl), [k] * len(pl), s=13, marker="s",
                     color="#e3e3e0", zorder=3, linewidths=0)
        if bad:
            ax_s.scatter(bad, [k] * len(bad), s=34, marker="s", color=C_DROWN,
                         zorder=5, linewidths=0)
        ax_s.text(len(DOMAINS) + 1.2, k,
                  f"{len(bad)} drown", va="center", ha="left", fontsize=7,
                  color=C_DROWN if bad else INK_MUTED)
    ax_s.set_yticks(range(len(rows)))
    ax_s.set_yticklabels([f"{METHOD_TICK[m]}\n{y}" for m, y in rows],
                         fontsize=7)
    ax_s.set_ylim(-0.6, len(rows) - 0.4)
    ax_s.set_xlim(0.5, len(DOMAINS) + 0.5)
    ax_s.set_xlabel(DOMAIN_AXIS_LABEL)
    _title(ax_s, 3, "domains that drown at initialisation")
    ax_s.grid(axis="x")
    ax_s.set_axisbelow(True)
    open_frame(ax_s)

    fig.legend(handles=[
        Line2D([], [], color=C_OLD, lw=2.8, label=METHOD_LABEL["old"]),
        Line2D([], [], color=C_NEW, lw=2.0, label=METHOD_LABEL["dunestart"]),
        Patch(facecolor=C_DROWN, label="drowns at initialisation"),
        Line2D([], [], color=NODATA, lw=8, label="outside the extraction"),
    ], loc="lower center", bbox_to_anchor=(0.5, -0.004), ncol=2,
        frameon=False, columnspacing=1.6, handlelength=2.4)

    med = "; ".join(f"{y} median {v:+.0f} m" for y, v in sorted(medians.items()))
    caption(fig, (
        "The two setback methods for NC-12, each period drawn on its own "
        "Barrier3D interiors. Domain 1 is at Cape Point in the south and "
        "domain 90 at Pea Island in the north; the shaded spans are the "
        "villages (Buxton, Avon, Tri-Village). The superseded method took the "
        "setback as the minimum road elevation minus the minimum dune "
        "elevation, independently per domain, referenced to the same-year "
        "digitised dune line; the current method measures the road landward "
        "from the dune start on each of a domain's 50 profiles and takes the "
        "median. Colour is the METHOD in this figure, not the vintage — the "
        "vintage is the panel. (a, b) both methods on one island each, the "
        "road drawn where roadway_manager.bulldoze puts it, road_start = "
        "int(setback / 10 m); topography "
        + ", ".join(f"{y} on {P.topo_label(y)}" for y in YEARS)
        + ", which are different islands rather than one island drawn twice. "
          "Interior elevation is in classes relative to mean high water; "
          "cells outside the extraction carry no data and are drawn grey. "
          "(c) the difference between the two, negative where the dune-start "
          f"method puts the road closer to the dune ({med}). (d) the domains "
          "where bulldoze's drown test fires at initialisation, one row per "
          "method and period; a red square is a roadway CASCADE stops "
          "managing."))

    save(fig, out_png)
    plt.close(fig)
    print(f"[out] {out_png}")


# =============================================================================
# FIGURE 2 -- both methods against the rasterized road
# =============================================================================

def figure_actual(per, crop_rows, placements, actual,
                  out_png: Path):
    fig = plt.figure(figsize=figsize("double", height=8.0))
    gs = fig.add_gridspec(4, 1, height_ratios=[1.25, 1.25, 1.0, 0.85],
                          hspace=0.30, left=0.135, right=0.870,
                          top=0.958, bottom=0.090)
    axes = [fig.add_subplot(gs[0])]
    axes += [fig.add_subplot(gs[i], sharex=axes[0]) for i in (1, 2, 3)]
    ax84, ax04, ax_e, ax_w = axes

    for i, (ax, year) in enumerate(((ax84, YEARS[0]), (ax04, YEARS[1]))):
        base_panel(ax, fig, per[year]["shown"], crop_rows,
                   f"the rasterized road in {year}, both methods over it", i)
        act = actual.get(year, {})
        if act:
            xs, lo, hi = [], [], []
            for d in sorted(act):
                a = act[d]
                xs += [d - 0.5, d + 0.5, np.nan]
                lo += [a["p10"], a["p10"], np.nan]
                hi += [a["p90"] + a["width"], a["p90"] + a["width"], np.nan]
            ax.fill_between(xs, lo, hi, color=C_ACTUAL, alpha=0.42, lw=0,
                            zorder=5)
        for name in METHOD_ORDER:
            if year in placements.get(name, {}):
                draw_method_line(ax, placements[name][year],
                                 METHOD_COLOUR[name], lw=2.0, zorder=7)

    # --- (C) error against the rasterized road ------------------------------
    # Both the per-DOMAIN error (line, against the domain's median road) and the
    # per-PROFILE spread of that error (band, against p10-p90 of the profiles).
    # The band is ported from the retired
    # 3-figures/island_wide/HAT_plot_road_placement_accuracy.py: a method can sit
    # on the median road and still miss most individual profiles, and only the
    # band shows that.
    ax_e.axhline(0, color=C_ACTUAL, lw=1.0, zorder=3)
    err_medians = {}
    for name in METHOD_ORDER:
        for year, style in zip(YEARS, [(0, ()), (0, (5, 1.6))]):
            if year not in placements.get(name, {}) or year not in actual:
                continue
            pl, act = placements[name][year], actual[year]
            common = sorted(set(pl) & set(act))
            err = [pl[d]["setback_m"] - act[d]["p50"] for d in common]
            err_medians[(name, year)] = float(np.median(err))
            if year == YEARS[0]:
                ax_e.fill_between(
                    common,
                    [pl[d]["setback_m"] - act[d]["p90"] for d in common],
                    [pl[d]["setback_m"] - act[d]["p10"] for d in common],
                    color=METHOD_COLOUR[name], alpha=0.16, lw=0, zorder=4)
            ax_e.plot(common, err, color=METHOD_COLOUR[name], lw=1.3, ls=style,
                      zorder=5, label=f"{METHOD_TICK[name]} {year}")
    ax_e.set_ylabel("method − rasterized\nroad (m)")
    ax_e.grid(axis="y")
    ax_e.set_axisbelow(True)
    town_bands(ax_e, label=False)
    open_frame(ax_e)
    ax_e.legend(loc="upper left", ncol=2, fontsize=7)
    _title(ax_e, 2, "distance from the road burnt on the grid")
    plt.setp(ax_e.get_xticklabels(), visible=False)

    # --- (D) what a scalar has to throw away --------------------------------
    spread_medians = {}
    for year, style in zip(YEARS, [(0, ()), (0, (5, 1.6))]):
        act = actual.get(year, {})
        if not act:
            continue
        xs = sorted(act)
        spread_medians[year] = float(np.median([act[d]["spread"] for d in xs]))
        ax_w.plot(xs, [act[d]["spread"] for d in xs], color=C_ACTUAL, lw=1.3,
                  ls=style, zorder=5, label=f"{year}")
    ax_w.set_ylabel("road position spread\nwithin a domain, p10–p90 (m)")
    ax_w.set_xlabel(DOMAIN_AXIS_LABEL)
    ax_w.set_xlim(0.5, len(DOMAINS) + 0.5)
    ax_w.set_ylim(bottom=0)
    ax_w.grid(axis="y")
    ax_w.set_axisbelow(True)
    town_bands(ax_w, label=False)
    open_frame(ax_w)
    ax_w.legend(loc="upper left", ncol=2, fontsize=7)
    _title(ax_w, 3, "how far the road moves within one domain")

    fig.legend(handles=[
        Patch(facecolor=C_ACTUAL, alpha=0.42,
              label="NC-12 as burnt on the domain grids (p10–p90 + width)"),
        Line2D([], [], color=C_OLD, lw=2.0, label=METHOD_LABEL["old"]),
        Line2D([], [], color=C_NEW, lw=2.0, label=METHOD_LABEL["dunestart"]),
        Line2D([], [], color=NODATA, lw=8, label="outside the extraction"),
    ], loc="lower center", bbox_to_anchor=(0.5, -0.004), ncol=2,
        frameon=False, columnspacing=1.6, handlelength=2.4)

    errs = "; ".join(f"{METHOD_TICK[m]} {y} median {v:+.0f} m"
                     for (m, y), v in sorted(err_medians.items()))
    spread = "; ".join(f"{y} median {v:.0f} m"
                       for y, v in sorted(spread_medians.items()))
    caption(fig, (
        "Both setback methods against NC-12 as it is actually rasterized onto "
        "the model grid. Domain 1 is at Cape Point in the south and domain 90 "
        "at Pea Island in the north; the shaded spans are the villages "
        "(Buxton, Avon, Tri-Village). The reference is the road burnt on the "
        "domain grids, p10–p90 of road_seaward_cell − interior_row0_cell "
        "across each domain's 50 profiles plus the measured road width, in "
        "the same frame both methods are drawn in. It is NOT an independent "
        "check on the dune-start method: that method's setback is the median "
        "of this band, so the two agree by construction up to int() "
        "truncation and the negative floor. It IS an independent check on the "
        "superseded method, which never saw this grid. Colour is the METHOD "
        "here, not the vintage. (a, b) the reference band with both methods "
        "over it, "
        + ", ".join(f"{y} on {P.topo_label(y)}" for y in YEARS)
        + ". Interior elevation is in classes relative to mean high water; "
          "cells outside the extraction carry no data and are drawn grey. "
          "(c) the distance from the domain's median rasterized road (line) "
          "and from its p10–p90 profiles (band, 1984 only): "
        + errs +
        ". (d) how far the road moves across a single domain's 50 profiles, "
        "the spread any scalar setback has to discard whichever method "
        "produced it: " + spread + "."))

    save(fig, out_png)
    plt.close(fig)
    print(f"[out] {out_png}")


# =============================================================================
# MAIN
# =============================================================================

def main() -> int:
    print("=" * 88)
    print("METHOD COMPARISON -- old vs dune-start, and both vs the real road")
    print("=" * 88)

    per, crop_rows, _max_rows = P.load_years(YEARS)
    for year in YEARS:
        print(f"  {year}: interiors from {P.topo_label(year)}")

    placements = load_placements(per)
    missing = [m for m in METHOD_ORDER if m not in placements]
    if missing:
        raise SystemExit(f"\n[stop] no setback files for: {missing}\n")
    actual = {y: load_actual(y) for y in YEARS}
    actual = {y: a for y, a in actual.items() if a}
    if not actual:
        print("  [warn] no profiles CSV -- figure 2 will have no reference band")

    for year in YEARS:
        if not all(year in placements[m] for m in METHOD_ORDER):
            continue
        o, n = placements["old"][year], placements["dunestart"][year]
        common = sorted(set(o) & set(n))
        diff = np.array([n[d]["setback_m"] - o[d]["setback_m"] for d in common])
        print(f"\n  {year}: dune-start - old  median {np.median(diff):+.0f} m | "
              f"p10 {np.percentile(diff, 10):+.0f} | "
              f"p90 {np.percentile(diff, 90):+.0f} | "
              f"max|{np.abs(diff).max():.0f}|")
        if year in actual:
            act = actual[year]
            for name in METHOD_ORDER:
                pl = placements[name][year]
                c = sorted(set(pl) & set(act))
                err = np.array([pl[d]["setback_m"] - act[d]["p50"] for d in c])
                print(f"    vs rasterized road, {METHOD_LABEL[name]:<18}: "
                      f"median {np.median(err):+6.0f} m | "
                      f"mean |err| {np.abs(err).mean():5.0f} m | "
                      f"max |err| {np.abs(err).max():5.0f} m")
            sp = np.array([act[d]["spread"] for d in act])
            print(f"    within-domain road spread (p10-p90): "
                  f"median {np.median(sp):.0f} m | max {sp.max():.0f} m")

    figure_methods(per, crop_rows, placements,
                   OUT_ROOT / "HAT_method_comparison_on_domains.png")
    figure_actual(per, crop_rows, placements, actual,
                  OUT_ROOT / "HAT_method_vs_actual_road.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
