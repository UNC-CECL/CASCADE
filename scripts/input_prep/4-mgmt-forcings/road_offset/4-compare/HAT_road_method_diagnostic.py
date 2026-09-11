# ==============================================================================
# HAT_road_method_diagnostic.py
#
# What changes when the road forcing moves from the legacy method to measuring
# against the extracted dune start -- and, crucially, WHY each domain changes.
#
# THE DECOMPOSITION
#   The legacy-to-new difference confounds two things. A control run with
#   STRAIGHTEN = False -- same code, same DEM dune crest, same road mask, same
#   median aggregation, only the frame changed -- separates them exactly:
#
#       new_straightened - legacy  =  (new_straightened - new_raw)   FRAME
#                                   + (new_raw          - legacy)    REFERENCE
#
#   FRAME      the obliquity correction. North-up clip boxes make NC-12 cross
#              each 500 m domain diagonally; straightening shears that out.
#   REFERENCE  everything else, and it is two effects that these files cannot
#              separate: the legacy method measures a digitized same-year
#              DUNE-LINE GEOJSON, this one measures the DEM DUNE CREST inside the
#              picked window, and those differ both in what feature they are and
#              in what year the feature dates from. Labelled honestly as one
#              component rather than split on an assumption.
#
#   Residual caveat: the two passes use two different pick files (a window picked
#   in one frame points at different cells in the other), so a little of FRAME is
#   really window choice. Second-order, not zero.
#
# WHAT IT SHOWS
#   A  setback per domain, three methods, per year
#   B  the two components
#   C  |total change| as a FRACTION OF ISLAND WIDTH -- a 50 m shift is severe on
#      a 150 m island and minor on a 600 m one, so metres alone hide which
#      domains matter
#   D  road elevation, legacy vs new: a near-null result, kept because it
#      documents that the elevation was never the problem
#
# USAGE
#     python HAT_road_method_diagnostic.py
# ==============================================================================

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Derived, not hardcoded: a machine-specific absolute path here made the script
# unrunnable for anyone else and silent about why. Two sibling scripts carried
# the same defect in a form that resolved to the filesystem root.
def _find_project_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit(f"cannot find data\\hatteras_init above {start}")


PROJECT_ROOT = _find_project_root(Path(__file__).resolve())
INIT_ROOT = PROJECT_ROOT / "data" / "hatteras_init"
# Topography version resolved from the extractor, not hardcoded -- it was
# "2009_v3" and silently survived the re-pick into 2009_v4. See hat_topo_version.py.
# parents[4] IS scripts/ -- hat_topo_version.py moved there 2026-08-20.
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
from hat_topo_version import (topo_dirs, array_name,  # noqa: E402
                             product_for_year)

# PER VINTAGE, not once (2026-08-26). This was a module-level topo_dirs() with
# no argument -- DEFAULT_PRODUCT -- and `island_width_m()` fed the SAME widths
# into both years of the table, so `total_change_frac_width` normalised the
# 1984 change by the 2004-start island. All 90 domains differ between the two
# products and 65 differ in interior SHAPE, so that denominator was wrong for
# every 1984 row. Same failure as the v3/v4 one above, in product form.
_TOPO_CACHE: dict[int, tuple] = {}


def topo_for_year(year: int):
    if year not in _TOPO_CACHE:
        _TOPO_CACHE[year] = topo_dirs(product_for_year(year))
    return _TOPO_CACHE[year]


def topo_label(year: int) -> str:
    return f"{product_for_year(year)}/{topo_for_year(year)[2]}"


# array_name() is the single definition of these filenames - the same one
# the extractor writes with. Nothing here spells a name.
ROADS_ROOT = INIT_ROOT / "4-mgmt-forcing" / "road_offset"
OFFSET_ROOT = ROADS_ROOT / "dunestart_offset"
LEGACY_SB_FMT = ROADS_ROOT / "old_method_offset" / "{year}" / "RoadSetback_{year}.csv"

# Road elevation is NOT per-year, and it stays that way -- but the reason is no
# longer "there is one 2009 DEM". There are two elevation products now, and in
# the road corridor they disagree by a median +0.222 m (2009-2014-1996 minus
# 2009-2014, 54 of 82 domains beyond 0.05 m). That difference is the
# uncorrected island-wide 1996-vs-2009 survey offset, not a roadbed, so it is
# kept OUT of the forcing: one set, sampled on the 2009-2014 baseline, read for
# both years in YEARS. The old per-year RoadElevation_<year>.csv pair does not
# exist -- writing two files implied a measured change in roadbed height
# between 1984 and 2004 that nothing supports.
# See data/.../road_elevation/RoadElevation_audit.md and the note beside
# HATTERAS_ROAD_ELEVATION_FILE in hatteras_site_config.py.
LEGACY_EL = (INIT_ROOT / "4-mgmt-forcing" / "road_elevation"
             / "RoadElevation.csv")

# 4-compare output lands in method_comparison/, NOT inside either method's
# folder -- this is a legacy-vs-dune-start comparison, so it belongs to neither.
# It wrote into dunestart_offset/ until 2026-08-28, which put a cross-method
# result under one of the methods it compares and contradicted the rule stated
# in this folder's README. Its sibling, HAT_method_comparison_figures.py, was
# always correct about that rule; only this one drifted.
OUT_ROOT = ROADS_ROOT / "method_comparison"
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
from hat_figure_style import (  # noqa: E402
    apply_style, figsize, save, caption, open_frame, town_bands, _title,
    DOMAIN_AXIS_LABEL, INK, INK_MUTED, C_1984 as C_EARLY, C_1997 as C_LATE,
)

apply_style()

OUT_PNG = OUT_ROOT / "HAT_road_method_diagnostic.png"
OUT_CSV = OUT_ROOT / "HAT_road_method_diagnostic.csv"

YEARS = [1984, 2004]
DOMAINS = list(range(1, 91))
CELL_SIZE_M = 10.0
SENTINEL_DAM = -0.3
BERM_MHW_M = 1.70 - 0.36

# Hue is the VINTAGE here and line style is the method, so the house pair
# applies: the earlier survey red, the later blue. Until 2026-09-10 this file
# had them the other way round (1984 blue, 2004 orange), which put it at odds
# with every other two-vintage figure in the project.
C_1984 = C_EARLY
C_2004 = C_LATE
INK_SECOND = INK
SURFACE = "white"

# Method is encoded by LINE STYLE, year by hue. That keeps the validated
# two-colour palette (normal-vision dE 33.6, worst CVD dE 26.5) and means no
# series is distinguished by colour alone.
# The prose that used to be burned onto the canvas: a title sentence, a
# three-line method paragraph and two statistics boxes. It goes to CAPTIONS.md
# beside the PNG, with the numbers filled from this run.
CAPTION_TMPL = (
    "Road forcing: what changes when the setback is measured from the extracted "
    "dune start instead of by the superseded method, and why. Island widths for "
    "the normalisation come from {topo_1984} and {topo_2004}. "
    "(a) Setback per domain under three methods, per year. "
    "(b) The difference split into its two parts. A control pass holds the DEM "
    "dune crest, the road mask, this code and the median aggregation fixed and "
    "changes only the frame, so the total splits exactly into FRAME (the "
    "obliquity correction: north-up clip boxes make NC-12 cross each 500 m "
    "domain diagonally) and REFERENCE (everything else). REFERENCE still mixes "
    "two effects these files cannot separate, a digitised same-year dune line "
    "against the DEM dune crest, and same-year against 2009, so it is reported "
    "as one component rather than split on an assumption. FRAME median "
    "{fr_med:+.0f} m, interquartile {fr_lo:+.0f} to {fr_hi:+.0f}, largest "
    "|{fr_max:.0f}|; REFERENCE median {rf_med:+.0f} m, largest |{rf_max:.0f}|. "
    "The change is almost entirely REFERENCE: straightening barely moves the "
    "per-domain scalar, because the median over 50 profiles already cancels a "
    "diagonal that is roughly symmetric about the domain centre. Obliquity "
    "still matters for the per-profile spread and for elevation sampling, just "
    "not for this median. "
    "(c) The same change as a fraction of island width, because a 50 m shift is "
    "severe on a 150 m island and minor on a 600 m one. "
    "(d) Road elevation under each method, a near-null result kept because it "
    "documents that elevation was never the problem. {elev}. "
    "CAVEAT, and it has grown: the two passes read two different pick files, so "
    "a little of FRAME is really window choice. Since the 2026-09-03 re-pick "
    "the production arm reads the current dune-topo version while the control "
    "arm still reads control-picks/HAT_dune_search_windows_1984-start_v1.json "
    "from 2026-08-27, so the two arms are now two VINTAGES as well as two "
    "frames, and panel (b)'s split carries that difference inside FRAME. "
    "Domain 1 is at Cape Point, domain 90 at Pea Island."
)


STYLE_NEW = (0, ())
STYLE_RAW = (0, (5, 1.6))
STYLE_LEGACY = (0, (1.4, 1.8))


def read_two_row(path: Path) -> dict[int, float]:
    if not path.is_file():
        return {}
    raw = np.loadtxt(path, delimiter=",")
    if raw.ndim != 2 or raw.shape[0] != 2:
        return {}
    return {int(k): float(v) for k, v in zip(raw[0], raw[1])}


def read_domains_csv(year: int, suffix: str = "") -> dict[int, dict]:
    path = OFFSET_ROOT / str(year) / f"RoadOffset_{year}_domains{suffix}.csv"
    if not path.is_file():
        return {}
    with open(path, newline="") as f:
        return {int(r["domain"]): r for r in csv.DictReader(f)}


def to_float(v) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return np.nan


def island_width_m(year: int) -> dict[int, float]:
    """Median land width per domain, from the interiors CASCADE will read.

    `year` is required: the two periods start from different extractions, so
    "the island width" is not one number per domain.
    """
    topo_dir = topo_for_year(year)[0]
    out = {}
    for d in DOMAINS:
        p = topo_dir / array_name("topography", d)
        if not p.is_file():
            continue
        topo = np.load(p)
        land = (topo > SENTINEL_DAM + 1e-6).sum(axis=0).astype(float)
        out[d] = float(np.median(land)) * CELL_SIZE_M
    return out


def build_table() -> list[dict]:
    rows = []
    for year in YEARS:
        widths = island_width_m(year)
        print(f"  {year}: island widths from {topo_label(year)}")
        new = read_domains_csv(year)
        raw = read_domains_csv(year, "_rawframe")
        leg_sb = read_two_row(Path(str(LEGACY_SB_FMT).format(year=year)))
        leg_el = read_two_row(LEGACY_EL)   # same file for every year in YEARS
        new_el = read_two_row(
            OFFSET_ROOT / str(year) / f"RoadElevation_{year}_dunestart.csv")

        for d in DOMAINS:
            if d not in new or int(new[d]["n_road_profiles"] or 0) == 0:
                continue
            # Honour the offset script's own span decision rather than
            # duplicating ROAD_SPAN here: domains it measured but excluded from
            # the model-facing files are not part of the forcing being compared.
            if "EXCLUDED_FROM_SPAN" in (new[d].get("flags") or ""):
                continue
            n_s = to_float(new[d].get("setback_dunestart_m"))
            n_r = (to_float(raw[d].get("setback_dunestart_m"))
                   if d in raw and int(raw[d]["n_road_profiles"] or 0) > 0
                   else np.nan)
            lg = leg_sb.get(d, np.nan)
            width = widths.get(d, np.nan)
            total = n_s - lg
            rows.append({
                "year": year, "domain": d,
                "section": new[d].get("section", ""),
                "setback_new_straight_m": n_s,
                "setback_new_raw_m": n_r,
                "setback_legacy_m": lg,
                "component_frame_m": n_s - n_r,
                "component_reference_m": n_r - lg,
                "total_change_m": total,
                "island_width_m": width,
                "total_change_frac_width": (abs(total) / width
                                            if width and np.isfinite(total)
                                            else np.nan),
                "elev_new_mhw": new_el.get(d, np.nan),
                "elev_legacy_mhw": leg_el.get(d, np.nan),
                "flags": new[d].get("flags", ""),
            })
    return rows


def series(rows, year, key):
    sel = [r for r in rows if r["year"] == year and np.isfinite(r[key])]
    return [r["domain"] for r in sel], [r[key] for r in sel]


def main() -> None:
    rows = build_table()
    if not rows:
        raise SystemExit("no diagnostic rows; run HAT_road_offset_from_dune_start "
                         "with CONTROL_UNSTRAIGHTENED = True first")

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    fig = plt.figure(figsize=figsize("double", height=9.0),
                     constrained_layout=True)
    gs = fig.add_gridspec(4, 1, height_ratios=[1.0, 1.0, 0.8, 0.8])
    ax_sb = fig.add_subplot(gs[0])
    ax_cp = fig.add_subplot(gs[1], sharex=ax_sb)
    ax_fr = fig.add_subplot(gs[2], sharex=ax_sb)
    ax_el = fig.add_subplot(gs[3], sharex=ax_sb)
    ax_el.set_xlim(0.5, len(DOMAINS) + 0.5)

    # ---------------- A: three methods ----------------
    ax_sb.axhline(0.0, color=INK, lw=0.8, zorder=3)
    for year, colour in ((1984, C_1984), (2004, C_2004)):
        for key, style, width in (("setback_new_straight_m", STYLE_NEW, 1.8),
                                  ("setback_new_raw_m", STYLE_RAW, 1.1),
                                  ("setback_legacy_m", STYLE_LEGACY, 1.1)):
            x, y = series(rows, year, key)
            ax_sb.plot(x, y, color=colour, lw=width, ls=style, zorder=6)
    ax_sb.set_ylabel("road setback\n(m landward of dune start)")
    _title(ax_sb, 0, "setback per domain, three methods")
    ax_sb.grid(axis="y")
    ax_sb.set_axisbelow(True)
    plt.setp(ax_sb.get_xticklabels(), visible=False)

    # ---------------- B: the two components ----------------
    ax_cp.axhline(0.0, color=INK, lw=0.8, zorder=3)
    for year, colour in ((1984, C_1984), (2004, C_2004)):
        x, y = series(rows, year, "component_reference_m")
        ax_cp.plot(x, y, color=colour, lw=1.8, ls=STYLE_NEW, zorder=6)
        x, y = series(rows, year, "component_frame_m")
        ax_cp.plot(x, y, color=colour, lw=1.1, ls=STYLE_RAW, zorder=7)

    fr = np.array([r["component_frame_m"] for r in rows
                   if np.isfinite(r["component_frame_m"])])
    rf = np.array([r["component_reference_m"] for r in rows
                   if np.isfinite(r["component_reference_m"])])
    ax_cp.set_ylabel("component of change\n(m)")
    _title(ax_cp, 1, "what the change is made of")
    ax_cp.grid(axis="y")
    ax_cp.set_axisbelow(True)
    plt.setp(ax_cp.get_xticklabels(), visible=False)

    # ---------------- C: normalized severity ----------------
    for year, colour, off in ((1984, C_1984, -0.2), (2004, C_2004, 0.2)):
        sel = [r for r in rows
               if r["year"] == year and np.isfinite(r["total_change_frac_width"])]
        ax_fr.bar([r["domain"] + off for r in sel],
                  [r["total_change_frac_width"] for r in sel],
                  width=0.4, color=colour, edgecolor="none", zorder=5)
    ax_fr.axhline(0.25, color=INK_MUTED, lw=0.8, ls=":", zorder=6)
    # Right-aligned: the tall bars are at GIS 11-15, so a left-hand label sits
    # on top of them.
    ax_fr.text(90.0, 0.26, "a quarter of the island's width",
               fontsize=7.5, color=INK_MUTED, va="bottom", ha="right")
    ax_fr.set_ylabel("|change| \u00f7 island width")
    _title(ax_fr, 2, "the same metres, against how wide the island is")
    ax_fr.grid(axis="y")
    ax_fr.set_axisbelow(True)
    plt.setp(ax_fr.get_xticklabels(), visible=False)

    # ---------------- D: elevation, the null result ----------------
    ax_el.axhline(BERM_MHW_M, color=INK_MUTED, lw=0.8, ls=":", zorder=3)
    ax_el.text(1.0, BERM_MHW_M + 0.03, f"berm {BERM_MHW_M:.2f} m MHW",
               fontsize=7.5, color=INK_MUTED, va="bottom")
    stats = []
    for year, colour in ((1984, C_1984), (2004, C_2004)):
        x, y = series(rows, year, "elev_new_mhw")
        ax_el.plot(x, y, color=colour, lw=1.8, ls=STYLE_NEW, zorder=6)
        x2, y2 = series(rows, year, "elev_legacy_mhw")
        ax_el.plot(x2, y2, color=colour, lw=1.1, ls=STYLE_LEGACY, zorder=5)
        pair = [(r["elev_legacy_mhw"], r["elev_new_mhw"]) for r in rows
                if r["year"] == year and np.isfinite(r["elev_legacy_mhw"])
                and np.isfinite(r["elev_new_mhw"])]
        a, b = np.asarray(pair, dtype=float).T
        stats.append(f"{year}: r = {np.corrcoef(a, b)[0, 1]:.3f}, "
                     f"median new minus legacy {np.median(b - a):+.02f} m, "
                     f"below the berm {int((b < BERM_MHW_M).sum())} of {len(b)} "
                     f"new against {int((a < BERM_MHW_M).sum())} of {len(a)} legacy")
    ax_el.set_ylabel("road elevation\n(m MHW)")
    _title(ax_el, 3, "road elevation, each method")
    ax_el.set_xlabel(DOMAIN_AXIS_LABEL)
    ax_el.grid(axis="y")
    ax_el.set_axisbelow(True)

    for ax in (ax_sb, ax_cp, ax_fr, ax_el):
        open_frame(ax)
        town_bands(ax, label=(ax is ax_sb))

    fig.legend(handles=[
        Line2D([], [], color=C_1984, lw=1.8, label="1984"),
        Line2D([], [], color=C_2004, lw=1.8, label="2004"),
        Line2D([], [], color=INK_MUTED, lw=1.8, ls=STYLE_NEW,
               label="measured from the dune start, straightened (production)"),
        Line2D([], [], color=INK_MUTED, lw=1.1, ls=STYLE_RAW,
               label="the same, unstraightened (control)"),
        Line2D([], [], color=INK_MUTED, lw=1.1, ls=STYLE_LEGACY,
               label="superseded method, from independent minima"),
    ], loc="outside lower center", ncol=3, frameon=False)

    caption(fig, CAPTION_TMPL.format(
        topo_1984=topo_label(1984), topo_2004=topo_label(2004),
        fr_med=np.median(fr), fr_lo=np.percentile(fr, 25),
        fr_hi=np.percentile(fr, 75), fr_max=np.abs(fr).max(),
        rf_med=np.median(rf), rf_max=np.abs(rf).max(),
        elev=" | ".join(stats)))
    save(fig, OUT_PNG, bbox_inches=None)
    plt.close(fig)

    print(f"[figure] {OUT_PNG}")
    print(f"[csv]    {OUT_CSV}  ({len(rows)} rows)")
    for year in YEARS:
        sel = [r for r in rows if r["year"] == year]
        tot = np.array([r["total_change_m"] for r in sel
                        if np.isfinite(r["total_change_m"])])
        frac = np.array([r["total_change_frac_width"] for r in sel
                         if np.isfinite(r["total_change_frac_width"])])
        big = [r["domain"] for r in sel
               if np.isfinite(r["total_change_frac_width"])
               and r["total_change_frac_width"] > 0.25]
        print(f"  {year}: total change median {np.median(tot):+.0f} m; "
              f"|change|/width median {np.median(frac):.0%}, "
              f"max {frac.max():.0%}; >25% of island width at {big}")


if __name__ == "__main__":
    main()
