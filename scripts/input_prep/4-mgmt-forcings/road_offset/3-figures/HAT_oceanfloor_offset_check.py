"""
How big the ocean-side move max(setback, 0) is, domain by domain, drawn where it fires.

    python scripts/input_prep/4-mgmt-forcings/road_offset/3-figures/HAT_oceanfloor_offset_check.py

Writes a CSV of the move per domain and zoomed figures of the two stretches. Details: scripts/input_prep/4-mgmt-forcings/road_offset/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-18
"""

from __future__ import annotations

import csv
import importlib.util
import math
import sys
import textwrap
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.patheffects import withStroke

# Shared code -- imported, never transcribed

PROJECT_ROOT = next(_p for _p in Path(__file__).resolve().parents
                    if (_p / "pyproject.toml").exists())
PLACEMENT = (PROJECT_ROOT / "scripts" / "input_prep" / "4-mgmt-forcings"
             / "road_offset" / "1-produce" / "HAT_road_placement_on_domains.py")


# HAT_road_placement_on_domains.py, loaded as a module
def load_placement():
    spec = importlib.util.spec_from_file_location("hat_placement", PLACEMENT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


P = load_placement()

INIT_ROOT = PROJECT_ROOT / "data" / "hatteras_init"
import sys as _tvsys
from pathlib import Path as _TVP
_tvsys.path.insert(0, str(next(_q for _q in _TVP(__file__).resolve().parents
                               if (_q / "pyproject.toml").exists()) / "scripts"))
from site_layer import hat_topo_version as _tv  # noqa: E402
# --- CONFIG ------------------------------------------------------------------
DUNESTART = _tv.ROAD_SETBACK_ROOT
OUT_DIR = DUNESTART / "modifications"
DOMAINS_FMT = "measured/{year}/RoadOffset_{year}_domains.csv"
PROFILES_FMT = "measured/{year}/RoadOffset_{year}_profiles.csv"

YEARS = (1984, 2004)
CELL = P.CELL_SIZE_M

C_SAND = "#e8dcc0"          # the strip seaward of interior row 0
C_RAW = "#8c0d24"           # measured, negative -- same crimson as the failures
C_MOVE = "#52514e"

CONTEXT_DOMAINS = 2         # unaffected neighbours drawn either side
CLUSTER_GAP = 3             # domains this close or closer share one zoom panel
HEADROOM_ABOVE_M = 170.0    # how far landward of row 0 the crop reaches
# -----------------------------------------------------------------------------


# Every in-span domain for a year, with both sides of the ocean-side move
def read_domains(year: int) -> list:
    path = DUNESTART / DOMAINS_FMT.format(year=year)
    if not path.is_file():
        return []
    out = []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            try:
                model = float(r["setback_model_m"])
                raw = float(r["setback_dunestart_m"])
                floored = float(r["setback_dunestart_floored_m"])
            except (KeyError, TypeError, ValueError):
                continue
            if not (math.isfinite(model) and math.isfinite(raw)):
                continue
            out.append(dict(domain=int(r["domain"]), section=r["section"],
                            raw=raw, floored=floored, model=model))
    return out


# Per-profile setbacks, keyed by domain, in profile order
def read_profiles(year: int) -> dict:
    path = DUNESTART / PROFILES_FMT.format(year=year)
    out = {}
    if not path.is_file():
        return out
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            try:
                out.setdefault(int(r["domain"]), []).append(
                    (int(r["profile"]), float(r["setback_m"])))
            except (KeyError, TypeError, ValueError):
                continue
    return {d: [v for _, v in sorted(vals)] for d, vals in out.items()}


# The affected domains and every number the report and figure need
def measure_move(year: int, interiors: dict) -> list:
    profiles = read_profiles(year)
    rows = []
    for rec in read_domains(year):
        raw = rec["raw"]
        if raw >= 0:
            continue
        d = rec["domain"]
        a = interiors.get(d)
        if a is None:
            continue
        n = a.shape[0]

        start_raw = int(raw / CELL)            # bulldoze's truncation
        end_raw = start_raw + P.ROAD_WIDTH_CELLS
        wrap_row = n + start_raw               # where numpy would really put it

        # The counterfactual: bulldoze's drown test on the wrapped road
        if abs(start_raw) < n:
            wrap_sea = P._wet_fraction(a[start_raw - 1, :])
            wrap_bay = P._wet_fraction(a[end_raw + 1, :])
        else:
            wrap_sea = wrap_bay = float("nan")

        # And at the destination, row 0, with the same test
        floor_sea = P._wet_fraction(a[0, :])
        floor_bay = P._wet_fraction(a[P.ROAD_WIDTH_CELLS + 1, :])

        ps = np.asarray(profiles.get(d, []), dtype=float)
        rows.append(dict(
            year=year, domain=d, section=rec["section"],
            raw_m=raw, floored_m=rec["floored"], model_m=rec["model"],
            delta_m=rec["floored"] - raw,
            delta_rows=int(rec["floored"] / CELL) - start_raw,
            island_m=n * CELL,
            pct_seaward=(float((ps < 0).mean() * 100) if ps.size
                         else float("nan")),
            n_profiles=int(ps.size),
            prof_min=(float(ps.min()) if ps.size else float("nan")),
            prof_max=(float(ps.max()) if ps.size else float("nan")),
            wrap_row=wrap_row, wrap_m=wrap_row * CELL,
            wrap_seaside=wrap_sea, wrap_bayside=wrap_bay,
            wrap_drowned=bool(wrap_sea > P.DROWN_PCT or wrap_bay > P.DROWN_PCT),
            floor_seaside=floor_sea, floor_bayside=floor_bay,
            floor_drowned=bool(floor_sea > P.DROWN_PCT
                               or floor_bay > P.DROWN_PCT),
            profiles=ps,
        ))
    return rows


# Every in-span domain's model-facing setback and profiles, for the year
def context_rows(year: int) -> dict:
    profiles = read_profiles(year)
    return {rec["domain"]: dict(model_m=rec["model"], raw_m=rec["raw"],
                                profiles=np.asarray(profiles.get(rec["domain"],
                                                                 []), float))
            for rec in read_domains(year)}


# Domains the move did NOT touch that still hold profiles seaward of row 0
def near_misses(year: int) -> list:
    profiles = read_profiles(year)
    out = []
    for rec in read_domains(year):
        if rec["raw"] < 0:
            continue
        ps = np.asarray(profiles.get(rec["domain"], []), dtype=float)
        if ps.size and (ps < 0).any():
            out.append((rec["domain"], float((ps < 0).mean() * 100)))
    return out


# Affected domains grouped into the stretches of island they sit on
def cluster(domains: list) -> list:
    out = []
    for d in sorted(set(domains)):
        if out and d - out[-1][-1] <= CLUSTER_GAP:
            out[-1].append(d)
        else:
            out.append([d])
    return out


# Figure

# One stretch of island, cropped to it, at a scale where 10 m is visible
def draw_zoom(ax, fig, canvas, rows_by_domain, context, lo, hi, floor_m,
              panel, show_cbar):
    ax.set_facecolor(P.WATER)
    c0, c1 = (lo - 1) * P.ALONG_COLS, hi * P.ALONG_COLS
    crop_rows = int((HEADROOM_ABOVE_M + CELL) / CELL)
    shown = canvas[:crop_rows, c0:c1]
    im = ax.imshow(np.ma.masked_invalid(shown), aspect="auto", origin="lower",
                   extent=[lo - 0.5, hi + 0.5, -CELL / 2,
                           crop_rows * CELL - CELL / 2],
                   cmap=P.LAND_CMAP, norm=Normalize(P.LAND_VMIN, P.LAND_VMAX),
                   interpolation="nearest")
    if show_cbar:
        cax = ax.inset_axes([1.015, 0.0, 0.018, 1.0])
        cb = fig.colorbar(im, cax=cax, extend="max")
        cb.set_label("interior elev (m MHW)", fontsize=8.5)
        cb.ax.tick_params(labelsize=8)
        cb.outline.set_edgecolor(P.INK_MUTED)

    # Seaward of interior row 0 the interior array holds nothing
    ax.axhspan(floor_m, -CELL / 2, color=C_SAND, lw=0, zorder=1)
    ax.axhline(-CELL / 2, color=P.INK_SECOND, lw=1.1, ls=(0, (4, 2)), zorder=5)

    halo = [withStroke(linewidth=4.0, foreground=P.SURFACE)]
    year = rows_by_domain[next(d for d in range(lo, hi + 1)
                               if d in rows_by_domain)]["year"]
    for d in range(lo, hi + 1):
        r = rows_by_domain.get(d)
        ctx = context.get(d)
        if r is None and ctx is None:
            continue                       # out of span -- nothing was measured

        # The model-facing band, all 20 m of it, drawn as a band rather than a line
        band_m = r["floored_m"] if r else ctx["model_m"]
        ax.add_patch(Rectangle((d - 0.5, band_m), 1.0,
                               P.ROAD_WIDTH_CELLS * CELL,
                               fc=P.C_YEAR[year], ec=P.SURFACE, lw=0.8,
                               alpha=0.95 if r else 0.34, zorder=8))

        # Every profile, at its own alongshore position, ON TOP of the band
        ps = (r or ctx)["profiles"]
        if ps.size:
            xs = d - 0.5 + (np.arange(ps.size) + 0.5) / ps.size
            neg = ps < 0
            ax.plot(xs[~neg], ps[~neg], lw=0, marker="o", ms=2.0,
                    mfc=P.INK_SECOND, mec="none", alpha=0.8, zorder=9)
            ax.plot(xs[neg], ps[neg], lw=0, marker="o", ms=2.4, mfc=C_RAW,
                    mec="none", alpha=0.9, zorder=9)
        if r is None:
            continue

        ax.plot([d - 0.5, d + 0.5], [r["raw_m"]] * 2, color=C_RAW, lw=3.0,
                solid_capstyle="butt", path_effects=halo, zorder=9)

        ax.annotate("", xy=(d, r["floored_m"]), xytext=(d, r["raw_m"]),
                    arrowprops=dict(arrowstyle="-|>", color=C_MOVE, lw=1.5,
                                    shrinkA=1.5, shrinkB=1.5,
                                    mutation_scale=13), zorder=10)
        ax.text(d + 0.07, (r["raw_m"] + r["floored_m"]) / 2,
                f"{r['delta_m']:.0f} m", ha="left", va="center", fontsize=9,
                color=C_MOVE, zorder=11,
                bbox=dict(fc=P.SURFACE, ec="none", alpha=0.82, pad=1.4))

    for (slo, shi), _ in P.SECTIONS:
        if lo <= shi + 0.5 <= hi:
            ax.axvline(shi + 0.5, color=P.SURFACE, lw=1.0, alpha=0.6, zorder=4)

    affected = [d for d in range(lo, hi + 1) if d in rows_by_domain]
    section = rows_by_domain[affected[0]]["section"]
    ax.set_title(f"({panel}) {section}, domains {affected[0]}–{affected[-1]}",
                 loc="left", fontsize=11, weight="semibold", pad=4)
    ax.set_xlim(lo - 0.5, hi + 0.5)
    ax.set_ylim(floor_m, HEADROOM_ABOVE_M)
    ax.set_xticks(range(lo, hi + 1))
    ax.set_xlabel("Barrier3D domain", fontsize=9)
    ax.set_ylabel("distance landward of interior row 0 (m)", fontsize=9)
    ax.grid(axis="y", color=P.INK_MUTED, alpha=0.18, lw=0.7)


# A per-domain table and a caption, in the register of a figure caption
def draw_caption(ax, rows, island_median, misses):
    ax.axis("off")
    head = (f"{'domain':>7}  {'measured':>10}  {'applied':>9}  "
            f"{'displacement':>14}  {'profiles seaward':>18}")
    unit = (f"{'':>7}  {'(m)':>10}  {'(m)':>9}  {'(m / rows)':>14}  "
            f"{'of row 0 (%)':>18}")
    lines = [head, unit, "-" * len(head)]
    for r in rows:
        lines.append(
            f"{r['domain']:>7}  {r['raw_m']:>10.0f}  {r['model_m']:>9.0f}  "
            f"{format(r['delta_m'], '.0f') + ' / ' + str(r['delta_rows']):>14}  "
            f"{r['pct_seaward']:>17.0f}%")
    ax.text(0.0, 1.0, "\n".join(lines), transform=ax.transAxes, va="top",
            ha="left", fontsize=9, family="monospace", color=P.INK_SECOND,
            linespacing=1.45)

    deltas = [r["delta_m"] for r in rows]
    caption = [
        f"Displacement is {min(deltas):.0f}–{max(deltas):.0f} m "
        f"({min(r['delta_rows'] for r in rows)}–"
        f"{max(r['delta_rows'] for r in rows)} grid rows; median "
        f"{np.median(deltas):.0f} m), against an island-wide median applied "
        f"setback of {island_median:.0f} m.",
        f"Left unconstrained, a negative setback is not rejected but indexes "
        f"the interior array from its landward margin, placing the roadway at "
        f"rows {min(r['wrap_row'] for r in rows)}–"
        f"{max(r['wrap_row'] for r in rows)} "
        f"({min(r['wrap_m'] for r in rows):.0f}–"
        f"{max(r['wrap_m'] for r in rows):.0f} m) in the back-barrier, where "
        f"{sum(r['wrap_drowned'] for r in rows)} of {len(rows)} fail the "
        f"initialisation drowning test; at row 0, "
        f"{sum(r['floor_drowned'] for r in rows)} of {len(rows)} fail.",
    ]
    if misses:
        caption.append(
            "The domain setback is the median of its cross-shore profiles, so "
            "the constraint is threshold-sensitive: "
            + "; ".join(f"domain {d} ({p:.0f}% of profiles seaward of row 0) is "
                        f"unconstrained" for d, p in misses)
            + f", while domain "
              f"{min(rows, key=lambda r: r['pct_seaward'])['domain']} "
              f"({min(r['pct_seaward'] for r in rows):.0f}%) is constrained.")
    # Wrapped here rather than with matplotlib's wrap=True
    ax.text(0.42, 1.0,
            "\n\n".join(textwrap.fill(c, 92) for c in caption),
            transform=ax.transAxes, va="top", ha="left", fontsize=9.2,
            color=P.INK_SECOND, linespacing=1.5)


# The check figure
def build_figure(rows: list, interiors: dict, island_median: float,
                 misses: list, panels_only: bool = False) -> Path:
    canvas, _ = P.build_canvas(interiors)
    by_domain = {r["domain"]: r for r in rows}
    groups = cluster([r["domain"] for r in rows])

    floor_m = min([r["raw_m"] for r in rows]
                  + [float(np.min(r["profiles"])) for r in rows
                     if r["profiles"].size]) - 25.0

    if panels_only:
        fig = plt.figure(figsize=(15.0, 6.0))
        gs = fig.add_gridspec(1, len(groups), wspace=0.15, left=0.058,
                              right=0.935, top=0.855, bottom=0.085)
    else:
        fig = plt.figure(figsize=(15.0, 9.0))
        gs = fig.add_gridspec(2, len(groups), height_ratios=[1.0, 0.34],
                              hspace=0.24, wspace=0.15, left=0.058,
                              right=0.935, top=0.800, bottom=0.055)
    years = sorted({r["year"] for r in rows})
    contexts = {y: context_rows(y) for y in years}
    for i, g in enumerate(groups):
        ax = fig.add_subplot(gs[0, i] if not panels_only else gs[i])
        draw_zoom(ax, fig, canvas, by_domain, contexts[by_domain[g[0]]["year"]],
                  max(1, g[0] - CONTEXT_DOMAINS),
                  min(len(P.DOMAINS), g[-1] + CONTEXT_DOMAINS),
                  floor_m, panel="abcdefgh"[i],
                  show_cbar=(i == len(groups) - 1))
    if not panels_only:
        draw_caption(fig.add_subplot(gs[1, :]), rows, island_median, misses)

    if not panels_only:
        fig.text(0.058, 0.982,
                 "Effect of the non-negative setback constraint on modelled NC-12 "
                 f"position, {', '.join(str(y) for y in years)}",
                 fontsize=14, va="top", weight="semibold")
        fig.text(0.058, 0.948,
                 f"Road setbacks were measured per cross-shore profile from the "
                 f"seaward edge of the Barrier3D interior (row 0) on the "
                 f"{', '.join(f'{y} on {P.topo_label(y)}' for y in YEARS)} "
                 f"— each vintage on its own DEM. "
                 f"{len(rows)} of {len(P.DOMAINS)} domains "
                 f"return a negative setback, i.e. the digitised roadway lies "
                 f"seaward of that boundary. CASCADE requires a\nnon-negative "
                 f"setback — `max(setback, 0)` — so these are displaced landward "
                 f"to row 0 before initialisation. Panels show the two affected "
                 f"alongshore reaches at the scale of the displacement; the "
                 f"whole-island view is in HAT_dunestart_stage0_raw.png.",
                 fontsize=9.2, color=P.INK_SECOND, va="top", linespacing=1.55)

    handles = [
        Line2D([], [], color=C_RAW, lw=3.0, label="measured setback"),
        Line2D([], [], color=P.C_YEAR[years[0]], lw=8,
               label="applied setback, 20 m road band"),
        Line2D([], [], color=P.C_YEAR[years[0]], lw=8, alpha=0.34,
               label="applied setback, unconstrained domain"),
        Line2D([], [], color=P.INK_SECOND, lw=0, marker="o", ms=5,
               label="individual profile (crimson: seaward of row 0)"),
        Line2D([], [], color=C_SAND, lw=8, label="seaward of interior row 0"),
    ]
    fig.legend(handles=handles, loc="upper left",
               bbox_to_anchor=(0.058, 0.985 if panels_only else 0.872),
               ncol=5, fontsize=9, framealpha=0.0, borderpad=0.4,
               columnspacing=1.6, handlelength=2.4)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_png = OUT_DIR / ("HAT_dunestart_oceanfloor_check_panels.png"
                         if panels_only
                         else "HAT_dunestart_oceanfloor_check.png")
    fig.savefig(out_png, dpi=150, bbox_inches="tight", facecolor=P.SURFACE)
    plt.close(fig)
    return out_png


# Report

CSV_FIELDS = ["year", "domain", "section", "raw_m", "floored_m", "model_m",
              "delta_m", "delta_rows", "pct_seaward", "n_profiles",
              "prof_min", "prof_max", "island_m", "wrap_row", "wrap_m",
              "wrap_seaside", "wrap_bayside", "wrap_drowned",
              "floor_seaside", "floor_bayside", "floor_drowned"]


# The per-domain move as CSV
def write_csv(rows: list) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / "HAT_dunestart_oceanfloor_check.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in CSV_FIELDS})
    return out


# Run: measure the move per vintage, write the CSV and figures
def main() -> int:
    # Per vintage: each year on its own island (2026-08-26)
    per, _crop_rows, _max_rows = P.load_years(YEARS)
    interiors_by_year = {y: per[y]["interiors"] for y in YEARS}
    rows, all_model, misses = [], [], []
    for year in YEARS:
        print(f"  {year}: interiors from {P.topo_label(year)}")
        rows += measure_move(year, interiors_by_year[year])
        all_model += [r["model"] for r in read_domains(year)]
        misses += near_misses(year)
    if not rows:
        print("no negative setbacks -- the ocean-side move fires nowhere")
        return 0
    island_median = float(np.median([v for v in all_model if v > 0]))

    print(f"{len(interiors_by_year[YEARS[0]])} interiors per vintage | "
          f"drown test imported from "
          f"{PLACEMENT.name}")
    print(f"\n{'=' * 98}")
    print("OCEAN-SIDE MOVE -- how different the model-facing offset is from "
          "the measurement")
    print("=" * 98)
    print(f"  {'year':>5} {'dom':>4} {'section':<26} {'measured':>9} "
          f"{'model':>7} {'move':>7} {'rows':>5} {'prof<0':>7} {'wrap':>10} "
          f"{'wrap drowns':>12}")
    for r in rows:
        print(f"  {r['year']:>5} {r['domain']:>4} {r['section'][:26]:<26} "
              f"{r['raw_m']:>8.0f}m {r['model_m']:>6.0f}m "
              f"{'+' + format(r['delta_m'], '.0f'):>6}m {r['delta_rows']:>5} "
              f"{r['pct_seaward']:>6.0f}% {r['wrap_m']:>8.0f}m "
              f"{'YES' if r['wrap_drowned'] else 'no':>12}")

    deltas = [r["delta_m"] for r in rows]
    print(f"\n  move        : median {np.median(deltas):.0f} m, "
          f"mean {np.mean(deltas):.0f} m, range "
          f"{min(deltas):.0f}-{max(deltas):.0f} m "
          f"({min(r['delta_rows'] for r in rows)}-"
          f"{max(r['delta_rows'] for r in rows)} model rows)")
    print(f"  in context  : island-wide median model setback is "
          f"{island_median:.0f} m, so the largest move is "
          f"{max(deltas) / island_median * 100:.0f}% of a typical setback")
    if misses:
        print(f"  sensitivity : the move fires on a median, so it is a vote. "
              f"Not floored, but holding negative profiles: "
              + ", ".join(f"domain {d} ({p:.0f}%)" for d, p in misses)
              + f" -- against {min(r['pct_seaward'] for r in rows):.0f}% at the "
                f"least-negative domain that WAS floored")
    print(f"  drown test  : {sum(r['wrap_drowned'] for r in rows)} of "
          f"{len(rows)} would drown at the wrapped row; "
          f"{sum(r['floor_drowned'] for r in rows)} of {len(rows)} drown at "
          f"row 0 where the move puts them")

    # The background strip is context for a per-domain measurement, not the measurement
    bg = interiors_by_year[YEARS[0]]
    print(f"\n  [out] {build_figure(rows, bg, island_median, misses)}")
    print(f"  [out] "
          f"{build_figure(rows, bg, island_median, misses, True)}")
    print(f"  [out] {write_csv(rows)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
