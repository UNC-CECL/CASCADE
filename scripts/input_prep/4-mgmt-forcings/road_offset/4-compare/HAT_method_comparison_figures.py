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

    blue    old method        (min road - min dune, digitised dune line)
    orange  dune-start method (per-profile, measured from interior row 0)

Same validated pair (OKLab normal-vision dE 33.6, worst CVD 26.5), reassigned.
Mixing the two conventions in one figure would be worse than reassigning them
in a figure that says so.

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
from matplotlib.colors import Normalize
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

ROADS_ROOT = P.ROADS_ROOT
OUT_ROOT = ROADS_ROOT
YEARS = P.YEARS
DOMAINS = P.DOMAINS
CELL_SIZE_M = P.CELL_SIZE_M
DROWN_PCT = P.DROWN_PCT
SURFACE, WATER = P.SURFACE, P.WATER
INK_MUTED, INK_SECOND = P.INK_MUTED, P.INK_SECOND
INK = "#0b0b0b"

# hue = METHOD here, not year. See the header.
C_OLD, C_NEW = "#2a78d6", "#eb6834"
C_ACTUAL = "#3b3a37"          # reference: neutral, never a category colour
C_DROWN = P.C_DROWN

METHOD_ORDER = ["old", "dunestart"]
METHOD_COLOUR = {"old": C_OLD, "dunestart": C_NEW}
METHOD_LABEL = {"old": "old method", "dunestart": "dune-start method"}

PROFILES_FMT = ("dunestart_offset/{year}/RoadOffset_{year}_profiles.csv")

plt.rcParams.update(P.plt.rcParams)


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

def base_panel(ax, fig, shown, crop_rows, title, label_sections):
    ax.set_facecolor(WATER)
    im = ax.imshow(np.ma.masked_invalid(shown), aspect="auto", origin="lower",
                   extent=[0.5, len(DOMAINS) + 0.5, -CELL_SIZE_M / 2,
                           crop_rows * CELL_SIZE_M - CELL_SIZE_M / 2],
                   cmap=P.LAND_CMAP, norm=Normalize(P.LAND_VMIN, P.LAND_VMAX),
                   interpolation="nearest")
    cax = ax.inset_axes([1.012, 0.0, 0.014, 1.0])
    cb = fig.colorbar(im, cax=cax, extend="max")
    cb.set_label("interior elev (m MHW)", fontsize=8.5)
    cb.ax.tick_params(labelsize=8)
    cb.outline.set_edgecolor(INK_MUTED)

    for (lo, hi), nm in P.SECTIONS:
        ax.axvline(hi + 0.5, color=SURFACE, lw=1.0, alpha=0.55, zorder=4)
        if label_sections:
            ax.text((lo + hi) / 2, crop_rows * CELL_SIZE_M * 0.97, nm,
                    ha="center", va="top", fontsize=8, color=INK_SECOND,
                    zorder=9,
                    bbox=dict(fc=SURFACE, ec="none", alpha=0.8, pad=1.6))

    ax.set_title(title, loc="left", fontsize=11.5, weight="semibold", pad=4)
    ax.set_ylabel("m landward of\ninterior row 0", fontsize=9)
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
    fig = plt.figure(figsize=(16.5, 13.4))
    gs = fig.add_gridspec(4, 1, height_ratios=[1.25, 1.25, 1.0, 0.62],
                          hspace=0.17, left=0.065, right=0.905,
                          top=0.872, bottom=0.055)
    axes = [fig.add_subplot(gs[0])]
    axes += [fig.add_subplot(gs[i], sharex=axes[0]) for i in (1, 2, 3)]
    ax84, ax04, ax_d, ax_s = axes

    for ax, year, lab in ((ax84, YEARS[0], True), (ax04, YEARS[1], False)):
        base_panel(ax, fig, per[year]["shown"], crop_rows,
                   f"{year} — both methods, on {P.topo_label(year)}", lab)
        for name in METHOD_ORDER:
            if year in placements.get(name, {}):
                draw_method_line(ax, placements[name][year],
                                 METHOD_COLOUR[name],
                                 lw=2.8 if name == "old" else 2.0,
                                 zorder=6 if name == "old" else 7)

    # --- (C) how far apart, per domain -------------------------------------
    ax_d.axhline(0, color=INK_SECOND, lw=1.1, zorder=3)
    for year, style in zip(YEARS, [(0, ()), (0, (5, 1.6))]):
        if not all(year in placements.get(m, {}) for m in METHOD_ORDER):
            continue
        o, n = (placements["old"][year], placements["dunestart"][year])
        common = sorted(set(o) & set(n))
        diff = [n[d]["setback_m"] - o[d]["setback_m"] for d in common]
        ax_d.plot(common, diff, color=INK, lw=1.7, ls=style, zorder=5,
                  label=f"{year}: dune-start − old  "
                        f"(median {np.median(diff):+.0f} m)")
    ax_d.set_ylabel("dune-start − old\nsetback (m)", fontsize=9)
    ax_d.grid(axis="y", color=INK_MUTED, alpha=0.22, lw=0.7)
    ax_d.set_axisbelow(True)
    ax_d.legend(loc="lower left", fontsize=8.5, ncol=2, framealpha=0.92)
    ax_d.set_title("(C)  Negative = the dune-start method puts the road CLOSER "
                   "to the dune", loc="left", fontsize=10)
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
                  f"{len(bad)} drown", va="center", ha="left", fontsize=8.5,
                  color=C_DROWN if bad else INK_MUTED)
    ax_s.set_yticks(range(len(rows)))
    ax_s.set_yticklabels([f"{METHOD_LABEL[m]}  {y}" for m, y in rows],
                         fontsize=8.5)
    ax_s.set_ylim(-0.6, len(rows) - 0.4)
    ax_s.set_xlim(0.5, len(DOMAINS) + 0.5)
    ax_s.set_xlabel("Barrier3D / GIS domain   "
                    "(1 = Cape Point / south  →  90 = Rodanthe / north)")
    ax_s.set_title("(D)  Drowns at initialisation — crimson square = CASCADE "
                   "stops managing this roadway", loc="left", fontsize=10)
    ax_s.grid(axis="x", color=INK_MUTED, alpha=0.18, lw=0.7)
    ax_s.set_axisbelow(True)

    fig.text(0.065, 0.985,
             "The two setback methods on each period's own Barrier3D "
             "interiors — where each one puts NC-12",
             fontsize=14, va="top", weight="semibold")
    fig.text(0.065, 0.958,
             "Topography: "
             + ", ".join(f"{y} on {P.topo_label(y)}" for y in YEARS)
             + ". Each panel compares two METHODS on ONE island; the two "
               "panels are different islands. HUE IS THE METHOD HERE, not "
               "the year — the year is the panel. Road drawn where bulldoze "
               "puts it: road_start = int(setback / 10 m).",
             fontsize=9, color=INK_SECOND, va="top", linespacing=1.5)
    fig.legend(handles=[
        Line2D([], [], color=C_OLD, lw=2.8, label="old method"),
        Line2D([], [], color=C_NEW, lw=2.0, label="dune-start method"),
        Patch(facecolor=C_DROWN, label="drowns at initialisation"),
        Line2D([], [], color=WATER, lw=8, label="off-island / sentinel water"),
    ], loc="upper left", bbox_to_anchor=(0.065, 0.933), ncol=4, fontsize=8.5,
        framealpha=0.0, borderpad=0.4, columnspacing=1.6, handlelength=2.6)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
    print(f"[out] {out_png}")


# =============================================================================
# FIGURE 2 -- both methods against the rasterized road
# =============================================================================

def figure_actual(per, crop_rows, placements, actual,
                  out_png: Path):
    fig = plt.figure(figsize=(16.5, 13.4))
    gs = fig.add_gridspec(4, 1, height_ratios=[1.25, 1.25, 1.0, 0.85],
                          hspace=0.17, left=0.065, right=0.905,
                          top=0.872, bottom=0.055)
    axes = [fig.add_subplot(gs[0])]
    axes += [fig.add_subplot(gs[i], sharex=axes[0]) for i in (1, 2, 3)]
    ax84, ax04, ax_e, ax_w = axes

    for ax, year, lab in ((ax84, YEARS[0], True), (ax04, YEARS[1], False)):
        base_panel(ax, fig, per[year]["shown"], crop_rows,
                   f"{year} — rasterized road, with both methods over it "
                   f"({P.topo_label(year)})", lab)
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
    ax_e.axhline(0, color=C_ACTUAL, lw=1.4, zorder=3)
    for name in METHOD_ORDER:
        for year, style in zip(YEARS, [(0, ()), (0, (5, 1.6))]):
            if year not in placements.get(name, {}) or year not in actual:
                continue
            pl, act = placements[name][year], actual[year]
            common = sorted(set(pl) & set(act))
            err = [pl[d]["setback_m"] - act[d]["p50"] for d in common]
            if year == YEARS[0]:
                ax_e.fill_between(
                    common,
                    [pl[d]["setback_m"] - act[d]["p90"] for d in common],
                    [pl[d]["setback_m"] - act[d]["p10"] for d in common],
                    color=METHOD_COLOUR[name], alpha=0.16, lw=0, zorder=4)
            ax_e.plot(common, err, color=METHOD_COLOUR[name], lw=1.7, ls=style,
                      zorder=5,
                      label=f"{METHOD_LABEL[name]} {year}  "
                            f"(median {np.median(err):+.0f} m)")
    ax_e.set_ylabel("method − rasterized\nroad (m)", fontsize=9)
    ax_e.grid(axis="y", color=INK_MUTED, alpha=0.22, lw=0.7)
    ax_e.set_axisbelow(True)
    ax_e.legend(loc="upper left", fontsize=8, ncol=2, framealpha=0.92)
    ax_e.set_title(f"(C)  Distance from the road actually burnt on the grid. "
                   f"Line = vs the domain's median road; band = vs its p10–p90 "
                   f"profiles ({YEARS[0]} only). The dune-start line is near "
                   f"zero BY CONSTRUCTION", loc="left", fontsize=10)
    plt.setp(ax_e.get_xticklabels(), visible=False)

    # --- (D) what a scalar has to throw away --------------------------------
    for year, style in zip(YEARS, [(0, ()), (0, (5, 1.6))]):
        act = actual.get(year, {})
        if not act:
            continue
        xs = sorted(act)
        ax_w.plot(xs, [act[d]["spread"] for d in xs], color=C_ACTUAL, lw=1.7,
                  ls=style, zorder=5,
                  label=f"{year}  (median {np.median([act[d]['spread'] for d in xs]):.0f} m)")
    ax_w.set_ylabel("road position spread\nwithin a domain, p10–p90 (m)",
                    fontsize=9)
    ax_w.set_xlabel("Barrier3D / GIS domain   "
                    "(1 = Cape Point / south  →  90 = Rodanthe / north)")
    ax_w.set_xlim(0.5, len(DOMAINS) + 0.5)
    ax_w.set_ylim(bottom=0)
    ax_w.grid(axis="y", color=INK_MUTED, alpha=0.22, lw=0.7)
    ax_w.set_axisbelow(True)
    ax_w.legend(loc="upper left", fontsize=8.5, ncol=2, framealpha=0.92)
    ax_w.set_title("(D)  How far the road moves across a single domain's 50 "
                   "profiles — the spread ANY scalar setback discards, "
                   "whichever method produced it", loc="left", fontsize=10)

    fig.text(0.065, 0.985,
             "Both methods against the road actually rasterized onto the model "
             "grid",
             fontsize=14, va="top", weight="semibold")
    fig.text(0.065, 0.958,
             "Reference band = NC-12 as burnt on the domain grids, p10–p90 "
             "across each domain's 50 profiles plus the measured road width, in "
             "the same frame both methods are drawn in.\n"
             "It is NOT an independent check on the dune-start method — that "
             "method's setback is the median of this band. It IS one on the old "
             "method, which never saw this grid.",
             fontsize=9, color=INK_SECOND, va="top", linespacing=1.5)
    fig.legend(handles=[
        Patch(facecolor=C_ACTUAL, alpha=0.42,
              label="rasterized NC-12 (p10–p90 + width)"),
        Line2D([], [], color=C_OLD, lw=2.0, label="old method"),
        Line2D([], [], color=C_NEW, lw=2.0, label="dune-start method"),
        Line2D([], [], color=WATER, lw=8, label="off-island / sentinel water"),
    ], loc="upper left", bbox_to_anchor=(0.065, 0.925), ncol=4, fontsize=8.5,
        framealpha=0.0, borderpad=0.4, columnspacing=1.6, handlelength=2.6)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150, bbox_inches="tight", facecolor=SURFACE)
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
