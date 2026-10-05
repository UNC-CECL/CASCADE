"""
When and where the beach was nourished: figures of the fill projects the hindcast fires.

    python scripts/input_prep/4-mgmt-forcings/beach_nourishment.py

Draws when and where each project fell, volume alongshore and a domain map,
from HATTERAS_NOURISHMENT_PROJECTS and the nourishment data. Details: scripts/input_prep/4-mgmt-forcings/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-18
"""

from __future__ import annotations

import datetime as dt
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = next(_p for _p in Path(__file__).resolve().parents
                    if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

from site_layer.hat_figure_style import (  # noqa: E402
    C, C_1984, C_1984_FILL, C_1997, C_1997_FILL, DOMAIN_AXIS_LABEL, INK,
    INK_MUTED, _halo, _title, apply_style, caption, figsize, open_frame,
    record_caption, save, town_bands,
)
from site_layer.hatteras_site_config import (  # noqa: E402
    HATTERAS_ANNOTATIONS, HATTERAS_DOMAINS, HATTERAS_NOURISHMENT_PROJECTS,
    HATTERAS_PERIODS,
)

# --- CONFIG ------------------------------------------------------------------
DATA_DIR = PROJECT_ROOT / "data" / "hatteras_init"
from site_layer.hat_topo_version import MGMT_ROOT as MGMT_DIR, MGMT_RECORD_XLSX as RECORD_XLSX  # noqa: E402
from site_layer.hat_observed_rates import DOMAIN_BOXES as DOMAIN_FILE  # noqa: E402
from site_layer.hat_map_layers import ISLAND_OUTLINE as OUTLINE_FILE  # noqa: E402
OUT_DIR = MGMT_DIR / "nourishment"

SPACING_M = HATTERAS_DOMAINS.domain_spacing_m
N_DOMAINS = HATTERAS_DOMAINS.num_real_domains
YEAR_LO, YEAR_HI = 1984, 2024

# Per-project maps
MAP_CRS = "EPSG:26918"          # UTM 18N, north up
MAP_PAD_DOMAINS = 3             # unfilled domains drawn either side of a footprint
MAP_ASPECT = 1.35               # map height / width
MAP_ZOOM = 15
GROIN_FILE = PROJECT_ROOT / "hard-structures" / "groin" / "HAT-groin-gis-analysis" / "gis_data" / "groins_hatteras.geojson"
TILE_CACHE = Path(__import__("tempfile").gettempdir()) / "hat_tile_cache"
# Maps are organised by place, so they colour by community (Okabe-Ito, clear of the year red/blue and the groin red)
# Okabe-Ito reddish purple, orange, yellow from north to south (2026-10-04): light and saturated so they
# read on aerial imagery, colour-blind safe, clear of the groin red and the locator teal. Drawn over a dark keyline
COMMUNITY_COLOUR = {"Rodanthe": ("#CC79A7", "#ebc4da"), "Avon": ("#E69F00", "#f5d48a"),
                    "Buxton": ("#F0E442", "#f8f1a6")}
KEYLINE = "#1a1a1a"
# Fills on the record but NOT in the model (they fall after the CoastSat data ends, 2026-01-13): maps only
from cascade_pipeline.nourishment import NourishmentProject  # noqa: E402
RECORD_ONLY_PROJECTS = (
    NourishmentProject(name="Avon 2026", year=2026, gis_domains=tuple(range(22, 27)),
                       volume_cubic_yards=375_000,
                       note="Pampas Drive (just south of Avon Pier) to Greenwood Place, ~1 mi; placed 2026-05-28 "
                            "to 2026-06-25. Greenwood Place geocodes ~190 m into GIS 22, the pier ~380 m into GIS 26"),
    NourishmentProject(name="Buxton 2026", year=2026, gis_domains=tuple(range(6, 17)),
                       volume_cubic_yards=2_000_000,
                       note="Haulover to the southernmost Buxton groin, ~2.9 mi; PLANNED volume, pumping from "
                            "2026-07-31, ~75% placed by 2026-09-16"),
)
PAPER_PAD_DOMAINS = 2           # unfilled domains either side on the paper figure
PAPER_ASPECT = 2.1              # paper map panel height / width

# Fill-year colours in time order: the earliest red, the latest blue, one between purple
_YEAR_PALETTE = [(C_1984, C_1984_FILL), (C["ACCENT"], C["ACCENT_FILL"]), (C_1997, C_1997_FILL)]
_FILL_YEARS = sorted({p.year for p in HATTERAS_NOURISHMENT_PROJECTS})
if len(_FILL_YEARS) > len(_YEAR_PALETTE):
    raise ValueError(f"{len(_FILL_YEARS)} fill years but {len(_YEAR_PALETTE)} colours; extend _YEAR_PALETTE")
_PICK = {1: [0], 2: [0, 2], 3: [0, 1, 2]}[len(_FILL_YEARS)]
YEAR_COLOUR = {y: _YEAR_PALETTE[i] for y, i in zip(_FILL_YEARS, _PICK)}
# -----------------------------------------------------------------------------


# One legend patch per fill year, named for its projects
def year_handles(lw, suffix=""):
    handles = []
    for y, (dark, light) in YEAR_COLOUR.items():
        names = [p.name.split()[0] for p in HATTERAS_NOURISHMENT_PROJECTS if p.year == y]
        handles.append(Patch(facecolor=light, edgecolor=dark, lw=lw,
                             label=f"{y}: {' + '.join(names)}{suffix}"))
    return handles


# Domains filled by more than one project, with the years
def refilled_domains():
    years = {}
    for p in HATTERAS_NOURISHMENT_PROJECTS:
        for g in p.gis_domains:
            years.setdefault(g, []).append(p.year)
    return {g: sorted(ys) for g, ys in years.items() if len(ys) > 1}


# The sentence naming each stretch filled more than once
def _refill_sentence():
    by_years = {}
    for g, ys in refilled_domains().items():
        by_years.setdefault(tuple(ys), []).append(g)
    times = {2: "twice", 3: "three times"}
    parts = [f"domains {min(gs)}–{max(gs)} are filled {times.get(len(ys), f'{len(ys)} times')} "
             f"({', '.join(map(str, ys))})"
             for ys, gs in by_years.items()]
    return "; ".join(parts)


# The two sources

# The site-config list, one row per project, with the volumes a run gets
def model_projects() -> pd.DataFrame:
    rows = []
    for p in HATTERAS_NOURISHMENT_PROJECTS:
        rows.append(dict(
            source="model input (hatteras_site_config)",
            name=p.name, year=p.year,
            first_gis=min(p.gis_domains), last_gis=max(p.gis_domains),
            n_domains=len(p.gis_domains),
            length_km=len(p.gis_domains) * SPACING_M / 1000,
            volume_cy=p.volume_cubic_yards,
            volume_m3=p.volume_m3_total,
            volume_m3_per_m=p.volume_m3_per_m(SPACING_M),
            in_modelled_reach=True, enabled=p.enabled, note=p.note,
        ))
    return pd.DataFrame(rows)


_VOL_RE = re.compile(r"([\d,]{5,})\s*cy")


# The spreadsheet as one row per year that carries a note or a flag
def record_projects() -> pd.DataFrame:
    raw = pd.read_excel(RECORD_XLSX, sheet_name="Nourishment_Timeline",
                        header=None, skiprows=3)
    rows = []
    for _, r in raw.iterrows():
        year = r.iloc[0]
        if not (isinstance(year, (int, float, np.integer, np.floating))
                and not pd.isna(year)):
            continue
        year = int(year)
        notes = " | ".join(str(v).strip() for v in r.iloc[1:3]
                           if isinstance(v, str) and v.strip())
        flags = pd.to_numeric(r.iloc[3:3 + N_DOMAINS], errors="coerce").to_numpy()
        gis = [i + 1 for i, v in enumerate(flags) if v == 1]
        if not gis and not notes:
            continue
        vols = [int(m.replace(",", "")) for m in _VOL_RE.findall(notes)]
        # One year can flag two stretches (2022: Buxton and Avon), so each run is its own row
        runs = _contiguous_runs(gis) or [[]]
        for run in runs:
            rows.append(dict(
                source="record (Hatteras_Management_Timelines.xlsx)",
                name=notes.split("(")[0].split("|")[0].strip() or f"{year} fill",
                year=year,
                first_gis=min(run) if run else np.nan,
                last_gis=max(run) if run else np.nan,
                n_domains=len(run),
                length_km=len(run) * SPACING_M / 1000 if run else np.nan,
                volume_cy=sum(vols) if vols else np.nan,
                volume_m3=np.nan, volume_m3_per_m=np.nan,
                in_modelled_reach=bool(run), enabled=np.nan, note=notes,
            ))
    return pd.DataFrame(rows)


# Consecutive domain numbers grouped into runs
def _contiguous_runs(numbers):
    runs = []
    for n in sorted(numbers):
        if runs and n == runs[-1][-1] + 1:
            runs[-1].append(n)
        else:
            runs.append([n])
    return runs


# (start, end) of every hindcast window in the site config
def period_windows():
    return sorted((int(k), int(v["end_year"])) for k, v in HATTERAS_PERIODS.items())


# Figure 1: when x where

# Each project by year and domain, against the hindcast windows
def fig_when_where(model: pd.DataFrame, record: pd.DataFrame) -> Path:
    fig, ax = plt.subplots(figsize=figsize("double", aspect=0.52),
                           constrained_layout=True)
    windows = period_windows()
    # The hindcast windows as brackets on the left, off-reach fills on the right
    x_lo, x_hi = -2.5 * len(windows) - 1.5, N_DOMAINS + 6.5
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(YEAR_LO - 0.8, YEAR_HI + 0.8)

    for y in range(YEAR_LO, YEAR_HI + 1, 5):
        ax.axhline(y, color=C["GRID"], lw=0.5, zorder=0)
    ax.axvline(0.5, color=INK_MUTED, lw=0.6, zorder=1)
    ax.axvline(N_DOMAINS + 0.5, color=INK_MUTED, lw=0.6, zorder=1)

    # Model extents: the light fill. Record flags: the darker inner bar.
    h_model, h_rec = 0.9, 0.42
    for _, p in model.iterrows():
        dark, light = YEAR_COLOUR[p.year]
        ax.add_patch(Rectangle((p.first_gis - 0.5, p.year - h_model / 2),
                               p.n_domains, h_model, facecolor=light,
                               edgecolor=dark, lw=0.7, zorder=3))
        ax.text((p.first_gis + p.last_gis) / 2, p.year + h_model / 2 + 0.25,
                f"{p['name'].split()[0]} {p.year}\n{p.first_gis}–{p.last_gis}",
                ha="center", va="bottom", fontsize=7.5, color=INK, zorder=6,
                path_effects=_halo(2.0))
    for _, r in record[record.in_modelled_reach].iterrows():
        dark, _ = YEAR_COLOUR.get(r.year, (INK, None))
        ax.add_patch(Rectangle((r.first_gis - 0.5, r.year - h_rec / 2),
                               r.n_domains, h_rec, facecolor=dark,
                               edgecolor="none", zorder=4))

    # Off-reach fills, one marker per year, north of GIS 90.
    off = record[~record.in_modelled_reach]
    x_off = N_DOMAINS + 3.5
    ax.scatter([x_off] * len(off), off.year, marker=">", s=22, facecolor="white",
               edgecolor=INK, lw=0.7, zorder=5)
    ax.text(x_off, YEAR_HI + 0.6, "north of\nthe reach", ha="center",
            va="bottom", fontsize=7, color=INK_MUTED, clip_on=False)

    # Hindcast windows.
    for k, (s, e) in enumerate(windows):
        x = -2.5 * (len(windows) - k) + 0.3
        ax.plot([x, x], [s, e], color=INK, lw=1.2, solid_capstyle="butt", zorder=3)
        for y in (s, e):
            ax.plot([x - 0.6, x + 0.6], [y, y], color=INK, lw=0.8, zorder=3)
        ax.text(x, (s + e) / 2, f"{s}–{e}", rotation=90, ha="center",
                va="center", fontsize=6.5, color=INK,
                bbox=dict(facecolor="white", edgecolor="none", pad=0.6), zorder=4)
    ax.text((x_lo + 0.5) / 2, YEAR_HI + 0.6, "hindcast\nwindows", ha="center",
            va="bottom", fontsize=7, color=INK_MUTED, clip_on=False)

    ax.set_xticks([1] + list(range(10, N_DOMAINS + 1, 10)))
    ax.set_yticks(range(YEAR_LO, YEAR_HI + 1, 5))
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    ax.set_ylabel("year placed")
    open_frame(ax)
    town_bands(ax, where="bottom", strip=0.045)

    handles = year_handles(0.7, ", model input extent") + [
        Patch(facecolor=INK, label="extent flagged in the management record"),
        Line2D([], [], marker=">", ls="none", markerfacecolor="white",
               markeredgecolor=INK, markersize=5,
               label="Pea Island / Oregon Inlet fill, north of GIS 90 (not modelled)"),
    ]
    fig.legend(handles=handles, loc="outside lower center", ncol=2, frameon=False)

    m = model.set_index("name")
    r = record[record.in_modelled_reach].set_index("year")
    n_off = int((~record.in_modelled_reach).sum())
    off_years = ", ".join(str(y) for y in off.year)
    # Model projects the record flags nothing for in their year
    unflagged = [f"{n} {int(row.year)}" for n, row in m.iterrows() if int(row.year) not in r.index]
    unflagged_text = (f" The record flags nothing for {', '.join(unflagged)}: the project is missing from "
                      "it (nourishment/datasets/README.md), so it has no inner bar." if unflagged else "")
    caption(fig, (
        f"Beach nourishment on the modelled reach, {YEAR_LO}–{YEAR_HI}, by year placed (vertical) "
        f"and GIS domain (horizontal, south to north, {SPACING_M:.0f} m each). Each filled bar is one "
        f"project as the hindcast receives it from hatteras_site_config.HATTERAS_NOURISHMENT_PROJECTS: "
        + "; ".join(f"{n} {int(row.year)}, domains {int(row.first_gis)}–{int(row.last_gis)} "
                    f"({row.n_domains} domains, {row.length_km:.1f} km, {row.volume_cy/1e6:.1f} M cy)"
                    for n, row in m.iterrows())
        + ". Colour is the year placed, red the earliest and blue the latest. The dark inner bar is the footprint flagged for "
          "the same project in Hatteras_Management_Timelines.xlsx (Nourishment_Timeline sheet): "
        + "; ".join(f"{int(y)} domains {int(row.first_gis)}–{int(row.last_gis)}" for y, row in r.iterrows())
        + ". Where the two differ the site-config extent was re-derived from the project description "
          "(Rodanthe: 2 mi north of the village, stopping short of the locked boundary domain 90; "
          "Avon: Due East Road to Askins Creek North Drive) and the record's flags are kept as drawn, "
          "not reconciled." + unflagged_text + " Hollow markers at the right are the {n_off} Pea Island / Oregon Inlet "
          f"navigation and emergency fills the record lists ({off_years}); they lie north of domain 90 "
          f"and no run sees them. Brackets at the left are the {len(windows)} hindcast windows; a run "
          "fires whatever falls inside its window, so "
        + _windows_sentence(windows, model.year.tolist())
        + ". Named bands along the foot are the villages."
    ).replace("{n_off}", str(n_off)))
    return save(fig, OUT_DIR / "nourishment_when_where", close=True)[0]


# The sentence saying which windows carry which fills, computed from the config so it cannot go stale
def _windows_sentence(windows, fill_years):
    def n_in(s, e):
        return sum(s <= y <= e for y in fill_years)
    by_count = {}
    for s, e in windows:
        by_count.setdefault(n_in(s, e), []).append(f"{s}–{e}")
    words = {0: "no fill", len(fill_years): f"all {len(fill_years)}"}
    parts = []
    for n in sorted(by_count):
        what = words.get(n, f"{n} of {len(fill_years)}")
        parts.append(" and ".join(by_count[n]) + (" carries " if len(by_count[n]) == 1 else " carry ") + what)
    return ", ".join(parts)


# Figure 2: volume alongshore

# Fill volume per domain, summed over the projects
def fig_volume_alongshore(model: pd.DataFrame) -> Path:
    fig, ax = plt.subplots(figsize=figsize("double", aspect=0.40),
                           constrained_layout=True)
    per_domain = np.zeros(N_DOMAINS + 1)
    refilled = refilled_domains()
    # Stacked in time order, so a refilled domain's later fill sits on its earlier one
    for p in sorted(HATTERAS_NOURISHMENT_PROJECTS, key=lambda q: q.year):
        v = p.volume_m3_per_m(SPACING_M)
        dark, light = YEAR_COLOUR[p.year]
        gis = np.array(p.gis_domains)
        base = per_domain[gis].copy()
        ax.bar(gis, [v] * len(gis), bottom=base, width=0.86, facecolor=light,
               edgecolor=dark, lw=0.6, zorder=3)
        per_domain[gis] += v
        # Above the bar, unless a later fill is stacked on it: then inside
        covered = any(refilled.get(int(g), [p.year])[-1] > p.year for g in gis)
        y, va = (base.max() + v / 2, "center") if covered else (base.max() + v + 12, "bottom")
        ax.text(gis.mean(), y, f"{p.name.split()[0]} {p.year}\n{v:.0f} m³/m",
                ha="center", va=va, fontsize=7.5, color=INK, zorder=6,
                path_effects=_halo(2.0))

    ax.set_xlim(0.5, N_DOMAINS + 0.5)
    ax.set_ylim(0, per_domain.max() * 1.32)
    ax.set_xticks([1] + list(range(10, N_DOMAINS + 1, 10)))
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    ax.set_ylabel("fill volume (m³ per m of shoreline)")
    ax.yaxis.grid(True, zorder=0)
    open_frame(ax)
    town_bands(ax, where="top", strip=0.07)

    # Fixed structures the footprints were measured from.
    for name, x in HATTERAS_ANNOTATIONS.groins.items():
        ax.axvline(x, color=INK, lw=0.8, ls=(0, (2, 2)), zorder=2)
        # South of the line: the footprints were measured north from it
        ax.text(x - 0.4, per_domain.max() * 0.92, name, rotation=90, ha="right",
                va="top", fontsize=6.5, color=INK_MUTED)
    for name, (x, _) in HATTERAS_ANNOTATIONS.piers.items():
        ax.plot([x], [0], marker="^", ms=5, color=INK, clip_on=False, zorder=7)
        ax.text(x + 0.4, per_domain.max() * 0.04, name, rotation=90, ha="left",
                va="bottom", fontsize=6.5, color=INK_MUTED, zorder=7)

    handles = year_handles(0.6) + [
        Line2D([], [], color=INK, lw=0.8, ls=(0, (2, 2)), label="groin"),
        Line2D([], [], marker="^", ls="none", color=INK, markersize=5, label="pier"),
    ]
    fig.legend(handles=handles, loc="outside lower center", ncol=len(handles), frameon=False)

    refill = _refill_sentence()
    refill_text = (f" Bars are stacked by year placed: {refill}, so the top of a stack is the "
                   f"{YEAR_LO}–{YEAR_HI} cumulative, {per_domain.max():.0f} m³/m at most."
                   if refill else f" No domain is filled twice, so the bars are also the "
                                  f"{YEAR_LO}–{YEAR_HI} cumulative.")
    dens = sorted(HATTERAS_NOURISHMENT_PROJECTS, key=lambda q: q.volume_m3_per_m(SPACING_M))
    lo, hi = dens[0], dens[-1]
    dens_text = (f" The densest fill, {hi.name.split()[0]} {hi.year}, put "
                 f"{hi.volume_m3_per_m(SPACING_M) / lo.volume_m3_per_m(SPACING_M):.1f} times as much "
                 f"sand per metre as the thinnest, {lo.name.split()[0]} {lo.year}.")

    parts = []
    for p in HATTERAS_NOURISHMENT_PROJECTS:
        parts.append(f"{p.name} ({p.year}): {p.volume_cubic_yards/1e6:.1f} M cy = "
                     f"{p.volume_m3_total/1e6:.2f} M m³ over {len(p.gis_domains)} domains, "
                     f"{p.volume_m3_per_m(SPACING_M):.0f} m³/m")
    caption(fig, (
        "Nourishment volume per metre of shoreline that each GIS domain receives in the hindcast, "
        "by year placed. A project's reported total (cubic yards, from the permitting record) is "
        f"converted to cubic metres and spread evenly over its domains of {SPACING_M:.0f} m: "
        + "; ".join(parts)
        + ". Even spreading is an assumption; real fill templates taper at their ends."
        + refill_text + dens_text
        + " The dotted line is the Buxton groin field, from which the Buxton footprint was measured "
          "north; triangles are the piers, and the Avon footprint was placed about the pier. Named "
          "bands along the top are the villages: the Buxton and Rodanthe footprints extend out of "
          "their villages into the road corridor, the Avon footprint stays inside its village."
    ))
    return save(fig, OUT_DIR / "nourishment_volume_alongshore", close=True)[0]


# Figure 3: the domain map

# The fills on a map of the rotated island
def fig_domain_map(model: pd.DataFrame) -> Path:
    import geopandas as gpd
    from shapely.affinity import rotate as shapely_rotate
    from shapely.ops import unary_union

    domains = gpd.read_file(DOMAIN_FILE)
    outline = gpd.read_file(OUTLINE_FILE).to_crs(domains.crs)
    origin = tuple(unary_union(domains.geometry.tolist()).centroid.coords[0])

    def to_strip(frame):
        # Rotate 90 deg clockwise so south is left, north right, ocean below
        return frame.set_geometry(frame.geometry.apply(
            lambda g: shapely_rotate(g, -90, origin=origin)), crs=frame.crs)

    strip = to_strip(domains)
    strip_outline = to_strip(outline)
    from shapely.geometry import box

    years_of = {}
    for p in sorted(HATTERAS_NOURISHMENT_PROJECTS, key=lambda q: q.year):
        for g in p.gis_domains:
            years_of.setdefault(g, []).append(p.year)
    year_of = years_of

    # A domain filled k times is cut into k cross-shore bands, the earliest fill on the ocean side
    pieces = []
    for _, row in strip[strip.domain_id.isin(years_of)].iterrows():
        ys = years_of[int(row.domain_id)]
        x0, y0, x1, y1 = row.geometry.bounds
        h = (y1 - y0) / len(ys)
        for k, y in enumerate(ys):
            pieces.append(dict(fill_year=y, geometry=row.geometry.intersection(
                box(x0, y0 + k * h, x1, y0 + (k + 1) * h))))
    bands = gpd.GeoDataFrame(pieces, geometry="geometry", crs=strip.crs)

    fig, ax = plt.subplots(figsize=figsize("double", aspect=0.36),
                           constrained_layout=True)
    minx, miny, maxx, maxy = strip.total_bounds
    strip_outline.plot(ax=ax, facecolor="0.93", edgecolor="none", zorder=0)
    strip.plot(ax=ax, facecolor="none", edgecolor="0.6", linewidth=0.35, zorder=1)
    for year, (dark, light) in YEAR_COLOUR.items():
        sub = bands[bands.fill_year == year]
        if not sub.empty:
            sub.plot(ax=ax, facecolor=light, edgecolor=dark, linewidth=0.9, zorder=3)

    dy = maxy - miny
    ax.set_xlim(minx - 0.01 * (maxx - minx), maxx + 0.01 * (maxx - minx))
    ax.set_ylim(miny - 0.55 * dy, maxy + 0.55 * dy)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)

    site_ends = set()
    for _, p in model.iterrows():
        site_ends |= {int(p.first_gis), int(p.last_gis)}
    # Every tenth domain, unless a footprint end sits within one domain of it (20 beside 21, 90 beside 89)
    label_domains = {d for d in range(10, N_DOMAINS + 1, 10)
                     if not any(abs(d - e) <= 1 for e in site_ends)} | {1} | site_ends
    for _, row in strip.iterrows():
        d = int(row["domain_id"])
        if d not in label_domains:
            continue
        pt = row.geometry.representative_point()
        is_site = d in year_of
        ax.annotate(str(d), xy=(pt.x, maxy + 0.04 * dy), ha="center", va="bottom",
                    fontsize=7, rotation=90, zorder=6,
                    color=INK if is_site else INK_MUTED,
                    fontweight="bold" if is_site else "normal")

    # One bracket per stretch: projects on the same domains share it
    spans = {}
    for _, p in model.sort_values("year").iterrows():
        spans.setdefault((int(p.first_gis), int(p.last_gis)), []).append(p)
    keys = sorted(spans)
    for a, b_ in zip(keys, keys[1:]):
        if b_[0] <= a[1]:
            raise ValueError(f"footprints {a} and {b_} overlap without matching; label them by hand")
    for (first, last), ps in spans.items():
        b = strip[strip.domain_id.between(first, last)].total_bounds
        ax.annotate("", xy=(b[0], maxy + 0.30 * dy), xytext=(b[2], maxy + 0.30 * dy),
                    arrowprops=dict(arrowstyle="|-|", linewidth=1.0, color=INK,
                                    mutation_scale=3), annotation_clip=False)
        ax.text((b[0] + b[2]) / 2, maxy + 0.34 * dy,
                f"{ps[0]['name'].split()[0]} {' + '.join(str(int(p.year)) for p in ps)}"
                f"\ndomains {first}–{last}",
                ha="center", va="bottom", fontsize=8, color=INK)

    ax.text(minx, miny - 0.16 * dy, "Cape Point (south)", ha="left", va="top",
            fontsize=8.5, color=INK_MUTED)
    ax.text(maxx, miny - 0.16 * dy,
            "Pea Island (north)\nOregon Inlet fills lie beyond this end, not modelled",
            ha="right", va="top", fontsize=8.5, color=INK_MUTED, linespacing=1.3)
    ax.text((minx + maxx) / 2, miny - 0.10 * dy, "Atlantic Ocean", ha="center",
            va="top", fontsize=8.5, color=INK_MUTED, style="italic")

    # A 5 km scale bar, and a north arrow along the rotated strip
    bx, by = minx + 0.02 * (maxx - minx), miny - 0.40 * dy
    ax.plot([bx, bx + 5000], [by, by], color=INK, lw=2.2, solid_capstyle="butt", zorder=12)
    ax.text(bx + 2500, by + 0.05 * dy, "5 km", ha="center", va="bottom", fontsize=8, color=INK)
    # Beside the scale bar, clear of the north-end labels.
    ax0, ax1 = bx + 7500, bx + 10000
    ax.annotate("", xy=(ax1, by), xytext=(ax0, by),
                arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.0, mutation_scale=11))
    ax.text((ax0 + ax1) / 2, by + 0.05 * dy, "N", ha="center", va="bottom",
            fontsize=8.5, fontweight="bold", color=INK)

    handles = year_handles(0.9) + [
        Patch(facecolor="white", edgecolor="0.6", lw=0.35, label="domain, never nourished"),
    ]
    fig.legend(handles=handles, loc="outside lower center", ncol=len(handles), frameon=False)
    refill = _refill_sentence()
    refill_text = (f" A domain filled more than once is split across the island into one band per "
                   f"fill, the earliest on the ocean side: {refill}." if refill else "")

    n_filled = len(year_of)
    caption(fig, (
        f"The {len(domains)} GIS domains in place, the island rotated 90° clockwise so that south is "
        "to the left, north to the right and the Atlantic at the bottom; rotation preserves distance, "
        f"so the scale bar holds. The {n_filled} domains that receive a nourishment project in the "
        "hindcast are filled by the year it was placed, with the site-config extents of Figure 1; "
        "grey is the island outline and white domains are never nourished." + refill_text + " Every tenth domain and "
        "the ends of each footprint are numbered along the top. The Pea Island / Oregon Inlet fills "
        "in the management record lie beyond the north end of the strip."
    ))
    return save(fig, OUT_DIR / "nourishment_domain_map", close=True)[0]


# The project maps and the summary table, in the paper style (2026-10-04): no text column, a panel
# label and one small corner box, lat/lon ticks, a locator, the house scale bar and north dart

# Short reported extent and source for each fill, for the table (full wording: reported_extent/reported_limits.csv)
REPORTED = {
    ("Rodanthe emergency fill", 2014): (
        "2.13 mi from 1.5 mi north of the Pea Island refuge border south into Mirlo Beach",
        "USACE public notice 2013, via Beachapedia",
        "https://beachapedia.org/State_of_the_Beach/State_Reports/NC/Beach_Fill"),
    ("Buxton beach nourishment", 2017): (
        "Haulover Day Use Area south to the groin at the old lighthouse site, 2.94 mi",
        "Outer Banks Voice, 2018-03-01",
        "https://outerbanksvoice.com/2018/03/01/delayed-buxton-beach-nourishment-project-is-finally-done/"),
    ("Buxton shore protection", 2022): (
        "Haulover Day Use Area to the lighthouse groin field, 2.9 mi",
        "Dare County bulletin, 2022",
        "https://content.govdelivery.com/accounts/NCDARECOUNTY/bulletins/3285a2f"),
    ("Avon shore protection", 2022): (
        "Due East Rd (3,000 ft north of Avon Pier) to the NPS/Avon boundary, 2.5 mi",
        "Dare County",
        "https://www.darenc.gov/government/beach-nourishment/avon-beach-nourishement"),
    ("Avon 2026", 2026): (
        "Just south of Avon Pier to Greenwood Place, about 1 mi",
        "Island Free Press, 2026",
        "https://islandfreepress.org/outer-banks-news/more-than-75-of-avon-nourishment-project-completed-as-work-moves-south/"),
    ("Buxton 2026", 2026): (
        "Haulover area to the southernmost Buxton groin, about 2.9 mi",
        "Island Free Press FAQ, 2026",
        "https://islandfreepress.org/blog/avon-and-buxton-beach-nourishment-faqs-2026-edition/"),
}
DATES_CSV = OUT_DIR / "datasets" / "nourishment_placement_dates.csv"
IMAGERY_NOTE = "Esri World Imagery, accessed 2026-10-04, not of the fill year"
MM = 1 / 25.4
PAPER_WIDTH_IN = 180 * MM       # journal double column
SINGLE_WIDTH_IN = 90 * MM       # about a single column
IMAGE_FADE = 0.15               # imagery blended this far toward white so the overlays read
PER_M = r"m$^{3}$ m$^{-1}$"     # mathtext, so the superscripts render in Arial


# Every fill, model ones first in time order, as (project, in_model)
def all_fills():
    every = [(q, True) for q in HATTERAS_NOURISHMENT_PROJECTS] + [(q, False) for q in RECORD_ONLY_PROJECTS]
    return sorted(every, key=lambda t: (t[0].year, not t[1], min(t[0].gis_domains)))


def _town(p):
    return p.name.split()[0]


# Placement dates per project name, from the dated record
def placement_dates():
    d = pd.read_csv(DATES_CSV).set_index("project")
    return {k: (r.start_date, r.end_date if isinstance(r.end_date, str) else "pending") for k, r in d.iterrows()}


# The map layers once: domains, their centres, groins and the island outline, all in MAP_CRS
def _layers():
    import geopandas as gpd
    dom = gpd.read_file(DOMAIN_FILE).to_crs(MAP_CRS).sort_values("domain_id").set_index("domain_id")
    return dict(dom=dom, cen=dom.geometry.centroid, groins=gpd.read_file(GROIN_FILE).to_crs(MAP_CRS),
                outline=gpd.read_file(OUTLINE_FILE).to_crs(MAP_CRS))


# A north-up window round a footprint with `pad` domains either side, at height/width `aspect`
def _window(dom, first, last, pad, aspect, east=0.0, min_h=0.0):
    lo, hi = max(1, first - pad), min(N_DOMAINS, last + pad)
    # Centred on the footprint itself, `pad` domains of margin above and below even at the ends of the reach
    x0, y0, x1, y1 = dom.loc[first:last].total_bounds
    y0, y1 = y0 - pad * SPACING_M, y1 + pad * SPACING_M
    x1 += east * (x1 - x0)
    cx_, cy_ = (x0 + x1) / 2, (y0 + y1) / 2
    hh = max(y1 - y0, (x1 - x0) * aspect, min_h) / 2 * 1.03
    return lo, hi, (cx_ - hh / aspect, cy_ - hh, cx_ + hh / aspect, cy_ + hh)


# Imagery for a window, faded toward white
def _imagery(ax, bx):
    import contextily as cx
    from rasterio.warp import transform_bounds
    cx.set_cache_dir(str(TILE_CACHE))
    w, s, e, n = transform_bounds(MAP_CRS, "EPSG:3857", *bx)
    img, ext = cx.bounds2img(w, s, e, n, zoom=MAP_ZOOM, source=cx.providers.Esri.WorldImagery, ll=False)
    img, ext = cx.warp_tiles(img, ext, t_crs=MAP_CRS)
    img = img.astype(float)
    img[..., :3] = img[..., :3] * (1 - IMAGE_FADE) + 255 * IMAGE_FADE
    ax.imshow(img.astype(np.uint8), extent=ext, zorder=0, interpolation="bilinear")


# Latitude ticks up the left edge and longitude ticks along the bottom, decimal degrees
def _latlon_ticks(ax, bx, fs, lat_step=0.02, lon_step=0.02):
    import pyproj
    to_ll = pyproj.Transformer.from_crs(MAP_CRS, "EPSG:4326", always_xy=True)
    to_utm = pyproj.Transformer.from_crs("EPSG:4326", MAP_CRS, always_xy=True)
    lon_l, lat_b = to_ll.transform(bx[0], bx[1])
    lon_r, lat_t = to_ll.transform(bx[2], bx[3])
    lats = np.arange(np.ceil(lat_b / lat_step) * lat_step, lat_t, lat_step)
    lons = np.arange(np.ceil(lon_l / lon_step) * lon_step, lon_r, lon_step)
    # Skip a tick that would sit in the frame corner, against the other axis's labels
    m_lat = 0.03 * (lat_t - lat_b)
    lats = lats[(lats > lat_b + m_lat) & (lats < lat_t - m_lat)]
    ys = [to_utm.transform(lon_l, la)[1] for la in lats]
    xs = [to_utm.transform(lo, lat_b)[0] for lo in lons]
    ax.set_yticks(ys)
    ax.set_yticklabels([f"{la:.2f}°N" for la in lats], rotation=90, va="center", fontsize=fs)
    ax.set_xticks(xs)
    ax.set_xticklabels([f"{-lo:.2f}°W" for lo in lons], fontsize=fs)
    ax.tick_params(length=2.0, width=0.5, pad=1.5)


# Panel label: the bold letter in its white box, the italic place and year beside it
def _panel_label(ax, i, place, fs, note=None):
    from site_layer.hat_figure_style import MAP_TEXT
    kw = {k: v for k, v in MAP_TEXT.items() if k != "fontsize"}
    kw["zorder"] = 20
    if i is None:
        t = ax.text(0.03, 0.975, place, transform=ax.transAxes, ha="left", va="top", fontsize=fs + 0.5, **kw)
    else:
        t0 = ax.text(0.03, 0.975, f"({chr(97 + i)})", transform=ax.transAxes, ha="left", va="top",
                     fontsize=fs + 0.5, fontweight="bold", **kw)
        t = ax.annotate(place, xy=(1, 0.5), xycoords=t0, xytext=(3, 0), textcoords="offset points", ha="left",
                        va="center", fontsize=fs + 0.5, **kw)
    if note:
        ax.annotate(note, xy=(0, 0), xycoords=t, xytext=(0, -2), textcoords="offset points", ha="left",
                    va="top", fontsize=fs - 1.5, fontstyle="italic", **kw)


# A plain scale bar: one black bar with a white keyline and a single centred label below
def _scale_bar(ax, length_m, at, fs):
    from site_layer.hat_figure_style import MAP_TEXT
    (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
    bx_, by_, h = x0 + at[0] * (x1 - x0), y0 + at[1] * (y1 - y0) + 0.025 * (y1 - y0), 0.008 * (y1 - y0)
    # Keyline the same weight as the panel border (spines are 0.5 pt)
    ax.add_patch(Rectangle((bx_, by_), length_m, h, facecolor=INK, edgecolor="white", lw=0.5, zorder=12))
    lab = f"{length_m / 1000:g} km" if length_m >= 1000 else f"{length_m:g} m"
    ax.text(bx_ + length_m / 2, by_ - 0.6 * h, lab, ha="center", va="top",
            **{**MAP_TEXT, "fontsize": fs, "zorder": 12})


# Village where it lies on the island: (GIS domain, fraction of the box width from its west edge)
VILLAGE_LABEL_AT = {"Rodanthe": (81, 0.55), "Buxton": (7.5, 0.35), "Avon": (29, 0.45)}


# Water bodies in spaced capitals and the village in italic, map-label style
def _geo_labels(ax, p, L, lo, hi, bx, fs, water=True, village=False):
    from site_layer.hat_figure_style import MAP_TEXT, place_label, water_label
    dom, w, h = L["dom"], bx[2] - bx[0], bx[3] - bx[1]
    kw = {**MAP_TEXT, "fontsize": fs - 1.0}
    # Beside the inset, not on it: up top when the inset is low, low when it is high
    # Single maps only: the paper panels are too narrow to hold them off the island (the locator orients those)
    if water:
        water_label(ax, bx[0] + 0.91 * w, bx[1] + 0.55 * h, "Atlantic Ocean", text_kw=kw, rotation=90)
        water_label(ax, bx[0] + 0.05 * w, bx[1] + 0.62 * h, "Pamlico Sound", text_kw=kw, rotation=90)
    # The village is named in the panel title, so it is not repeated on the map (Hannah, 2026-10-04)
    if not village:
        return
    town = _town(p)
    d, f = VILLAGE_LABEL_AT[town]
    if town == "Avon" and not lo <= d <= hi:
        d = hi - 0.5
    if lo <= d <= hi:
        b = dom.loc[int(d)].geometry.bounds
        y = b[1] + (d - int(d) + 0.5 if d != int(d) else 0.5) * (b[3] - b[1])
        place_label(ax, b[0] + f * (b[2] - b[0]), y, town, text_kw={**MAP_TEXT, "fontsize": fs - 1.0})


# One fill on its imagery: domains, footprint, end numbers, groins and piers, scale bar and north dart
def _draw_fill(ax, p, in_model, L, pad, aspect, i, fs, corner=(0.97, 0.975, "right", "top"), oneline=False,
               east=0.0, north_top=False, north=True, title_outside=False, scale_at=(0.07, 0.065),
               years=None, scale_m=500, north_left=False, min_h=0.0):
    from shapely.ops import unary_union
    from site_layer.hat_figure_style import MAP_TEXT, north_dart, scale_bar_km
    dom, cen = L["dom"], L["cen"]
    first, last = min(p.gis_domains), max(p.gis_domains)
    dark, light = COMMUNITY_COLOUR[_town(p)]
    lo, hi, bx = _window(dom, first, last, pad, aspect, east, min_h)
    _imagery(ax, bx)
    # Only the footprint's own domains are drawn; the unfilled ones either side are left off
    dom.loc[first:last].plot(ax=ax, facecolor="none", edgecolor="white", lw=0.3, alpha=0.6, zorder=2)
    # Merged into one outline for the whole section; a wider buffer closes the small gaps between domain boxes
    fp = unary_union([g.buffer(25) for g in dom.loc[first:last].geometry]).buffer(-25)
    import geopandas as gpd
    gpd.GeoSeries([fp], crs=MAP_CRS).plot(ax=ax, facecolor=light, edgecolor="none", alpha=0.15, zorder=3)
    edge = gpd.GeoSeries([fp.exterior if fp.geom_type == "Polygon" else fp.boundary], crs=MAP_CRS)
    ls = "-" if in_model else (0, (3, 1.5))
    # A dark keyline under the colour so the outline reads on pale beach and dark water alike
    edge.plot(ax=ax, color=KEYLINE, lw=2.4, zorder=4, linestyle=ls)
    edge.plot(ax=ax, color=dark, lw=1.2, zorder=4.1, linestyle=ls)
    for d in (first, last):
        b = dom.loc[d].geometry.bounds
        ax.text(b[0] + 70, cen.loc[d].y, str(d), ha="left", va="center", **{**MAP_TEXT, "fontsize": fs - 1.0})
    hit = L["groins"][L["groins"].intersects(dom.loc[lo:hi].union_all())]
    ax._has_groins = bool(len(hit))
    ax._has_pier = any(lo <= d < hi for d, _ in HATTERAS_ANNOTATIONS.piers.values())
    if len(hit):
        hit.plot(ax=ax, color="white", lw=2.6, zorder=7)
        hit.plot(ax=ax, color=C["GROIN"], lw=1.3, zorder=7)
    for name, (d, frac) in HATTERAS_ANNOTATIONS.piers.items():
        if lo <= d < hi:
            a, b2 = cen.loc[d], cen.loc[d + 1]
            py = a.y + (frac - 0.5) * (b2.y - a.y)
            px = dom.loc[d].geometry.bounds[2] - 200
            # Drawn like the groins: a short shore-normal line, white with a dark keyline, ~220 m seaward
            t = np.array([b2.x - a.x, b2.y - a.y]); t /= np.hypot(*t)
            nrm = np.array([t[1], -t[0]]) if t[1] > 0 else np.array([-t[1], t[0]])
            p0, p1 = np.array([px, py]) - 40 * nrm, np.array([px, py]) + 220 * nrm
            ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=INK, lw=3.0, solid_capstyle="butt", zorder=7)
            ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color="white", lw=1.4, solid_capstyle="butt", zorder=7)
    ax.set_xlim(bx[0], bx[2])
    ax.set_ylim(bx[1], bx[3])
    ax.set_aspect("equal")
    for sp in ax.spines.values():
        sp.set_linewidth(0.5)
        sp.set_color(INK)
    _latlon_ticks(ax, bx, fs - 1.5)
    if title_outside:
        # Panel title above the frame, centred: the letter and the place
        ax.set_title(f"({chr(97 + i)}) {_town(p)}, {years or p.year}", fontsize=fs + 0.5, color=INK, pad=3)
    else:
        _panel_label(ax, i, f"{_town(p)}, {p.year}", fs, note=None if in_model else "not in the model")
    _geo_labels(ax, p, L, lo, hi, bx, fs, water=i is None)
    # A short bar, placed on open water by the caller
    # "mirror": the lower-left arrow-and-bar layout flipped to the lower right, same spacing from the edge
    if north_left == "mirror":
        scale_at = (1 - 0.17 - scale_m / (bx[2] - bx[0]), 0.03)
    _scale_bar(ax, scale_m, scale_at, fs - 1.0)
    if north:
        w_ = bx[2] - bx[0]
        nx = (bx[2] - 0.08 * w_ if north_left == "mirror" else bx[0] + 0.08 * w_ if north_left
              else bx[2] - 0.13 * w_)
        north_dart(ax, (nx, bx[1] + (0.86 if north_top else 0.07 if north_left else 0.10) * (bx[3] - bx[1])),
                   arrow_m=0.055 * (bx[3] - bx[1]), text_kw={**MAP_TEXT, "fontsize": fs - 1.0})
    return bx


# The locator: the island outline with each panel's window boxed and lettered
# Every village on the reach, nourished or not, at its GIS domain (north to south)
LOCATOR_VILLAGES = {"Rodanthe": 80, "Waves": 74, "Salvo": 69, "Avon": 26, "Buxton": 7.5}


def _locator(ax, L, windows, fs, label_side="right", fill=False, villages=False):
    from site_layer.hat_figure_style import place_label
    # Zoomed to the 90 model domains, not the whole outline: land grey on pale water, the reach outlined
    ax.set_facecolor("#eef3f6")
    L["outline"].plot(ax=ax, facecolor="0.82", edgecolor="0.45", lw=0.35, zorder=1)
    reach = L["dom"].union_all()
    import geopandas as gpd
    gpd.GeoSeries([reach.envelope.buffer(0)], crs=MAP_CRS).boundary.plot(ax=ax, color="none", lw=0)
    ob = L["dom"].total_bounds
    merged = {}
    for letter, bx in windows:
        merged.setdefault(tuple(np.round(bx, -1)), (bx, []))[1].append(letter)
    for bx, letters in merged.values():
        letter = ", ".join(x for x in letters if x)
        ax.add_patch(Rectangle((bx[0], bx[1]), bx[2] - bx[0], bx[3] - bx[1], facecolor="none",
                               edgecolor=C["LOCATOR"], lw=0.9, zorder=3))
        if letter:
            xl = bx[2] + 1500 if label_side == "right" else bx[0] - 1500
            ax.text(xl, (bx[1] + bx[3]) / 2, letter, ha="left" if label_side == "right" else "right",
                    va="center", fontsize=fs, color=C["LOCATOR"], fontweight="bold", zorder=4)
    # Zoomed to the panels it locates (plus a little), not the whole 90-domain reach
    wb = np.array([b for _, b in windows])
    # Down past the cape point when villages are labelled, so the cape is not cut off
    # With villages: down past the cape point, and open water above so the island sits lower in the frame
    y_lo, y_hi = wb[:, 1].min() - (3500 if villages else 1500), wb[:, 3].max() + (3000 if villages else 1500)
    if villages:
        y_lo, y_hi = y_lo - 2000, y_hi - 2000          # view nudged south so the island sits a little higher
    d = L["dom"].geometry.bounds
    near = d[(d.maxy > y_lo) & (d.miny < y_hi)]
    if villages:
        # Centred on the boxes; equal aspect then sets the width from the panel's shape
        xc = (wb[:, 0].min() + wb[:, 2].max()) / 2 - 1200   # nudged west so the island sits a little right
        ax.set_xlim(xc - 2000, xc + 2000)
    else:
        ax.set_xlim(min(near.minx.min(), wb[:, 0].min()) - 500, max(near.maxx.max(), wb[:, 2].max()) + 3500)
    ax.set_ylim(y_lo, y_hi)
    if label_side == "right" and not any(x for x, _ in windows):
        place_label(ax, ob[0] + 0.35 * (ob[2] - ob[0]), ob[1] + 0.55 * (ob[3] - ob[1]), "Hatteras\nIsland",
                    text_kw=dict(color=INK, fontsize=fs - 1.0, zorder=5))
    # Small dots and names on the sound side for every village, so the boxes sit among their neighbours
    if villages:
        from shapely.geometry import LineString
        cen = L["cen"]
        land = L["outline"].union_all()
        for name, g in LOCATOR_VILLAGES.items():
            pt = cen.loc[int(g)]
            b = L["dom"].loc[int(g)].geometry.bounds
            if not y_lo <= pt.y <= y_hi:
                continue
            # The dot goes on the island: the middle of the widest land crossing at the village's northing
            cut = LineString([(b[0] - 3000, pt.y), (b[2] + 500, pt.y)]).intersection(land)
            parts = list(getattr(cut, "geoms", [cut]))
            seg = max(parts, key=lambda q: q.length) if parts and not cut.is_empty else None
            xd = seg.centroid.x if seg is not None else (b[0] + b[2]) / 2
            ax.plot(xd, pt.y, marker="o", ms=1.8, color=INK, zorder=5)
            # Buxton's name drops below its dot, onto the cape
            dy = -1900 if name == "Buxton" else 0
            # Buxton's name sits centred under its dot, on the cape, clear of the panel edge
            tx, ha = (xd, "center") if dy else (b[0] - 600, "right")
            ax.text(tx, pt.y + dy, name, ha=ha, va="center", fontsize=fs - 1.0, fontstyle="italic",
                    color=INK, zorder=5)
    # fill=True keeps the inset box exactly where it was placed and widens the view instead of shrinking the box
    ax.set_aspect("equal", adjustable="datalim" if fill else "box")
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_linewidth(0.5)
        sp.set_color(INK)


def _legend_handles(towns, record_only=False, groins=True, pier=True):
    import matplotlib.patheffects as pe
    h = [Patch(facecolor=COMMUNITY_COLOUR[t][1], edgecolor=COMMUNITY_COLOUR[t][0], lw=1.2, label=t,
               path_effects=[pe.Stroke(linewidth=2.4, foreground=KEYLINE), pe.Normal()])
         for t in towns]
    if record_only:
        h.append(Line2D([], [], color=INK, lw=0.9, ls=(0, (3, 1.5)), label="Not in the model"))
    if groins:
        h.append(Line2D([], [], color=C["GROIN"], lw=1.3, label="Buxton groins"))
    if pier:
        h.append(Line2D([], [], color="white", lw=1.4, label="Pier",
                        path_effects=[__import__("matplotlib.patheffects", fromlist=["Stroke"]).Stroke(
                            linewidth=3.0, foreground=INK), __import__("matplotlib.patheffects",
                            fromlist=["Normal"]).Normal()]))
    h.append(Patch(facecolor="none", edgecolor="0.6", lw=0.5, label="Model domain"))
    return h


def _caption_row(p, in_model, dates):
    per_m = p.volume_m3_per_m(SPACING_M)
    s, e = dates.get(p.name, ("?", "?"))
    return (f"{p.name}, {p.year}: GIS {min(p.gis_domains)}–{max(p.gis_domains)} "
            f"({len(p.gis_domains) * SPACING_M / 1000:.1f} km), {p.volume_cubic_yards / 1e6:.2f} × 10⁶ yd³ "
            f"({p.volume_m3_total / 1e6:.2f} × 10⁶ m³, {per_m:.0f} m³ m⁻¹), placed {s} to {e}"
            + ("" if in_model else ", on the record but not in the model"))


# One map per fill, the supplement versions: same style as the paper panels, locator in the upper right
def fig_project_maps(model: pd.DataFrame) -> list[Path]:
    L = _layers()
    dates = placement_dates()
    outs = []
    aspect = 1.55
    # One map scale for every single map, so the north arrow and the 1 km bar are the same size on each
    min_h = max(L["dom"].loc[min(q.gis_domains):max(q.gis_domains)].total_bounds[3]
                - L["dom"].loc[min(q.gis_domains):max(q.gis_domains)].total_bounds[1]
                for q, _ in all_fills()) + 2 * MAP_PAD_DOMAINS * SPACING_M
    for p, in_model in all_fills():
        fig = plt.figure(figsize=(SINGLE_WIDTH_IN, SINGLE_WIDTH_IN * aspect * 0.86))
        ax = fig.add_axes([0.13, 0.11, 0.84, 0.86])
        # Info box under the place label; the locator takes the upper right, over the ocean, clear of the island
        bx = _draw_fill(ax, p, in_model, L, MAP_PAD_DOMAINS, aspect, None, fs=8.0,
                        scale_at=(0.15, 0.03), scale_m=1000, min_h=min_h,
                        north_left=True)
        # Where the beach runs up the upper right the inset sits lower right and the north arrow goes up top
        # One inset size for every map (same figure and axes size, so the same box in inches)
        # Tucked into the corner, a hair off both edges
        # No locator on the single maps: they show the region only (Hannah, 2026-10-04); the paper figure has one
        hs = _legend_handles([_town(p)], record_only=not in_model, groins=ax._has_groins, pier=ax._has_pier)
        # Legend inside the map, upper right over the ocean, in the house map-legend box
        from site_layer.hat_figure_style import MAP_LEGEND
        # Avon and Buxton: legend lower right, clear of the island at the top of the frame
        ax.legend(handles=hs, loc="lower right" if _town(p) in ("Avon", "Buxton") else "upper right",
                  **{**MAP_LEGEND, "fontsize": 7, "borderaxespad": 0.4})
        stem = f"nourishment_map_{p.year}_{_town(p).lower()}"
        out = save(fig, OUT_DIR / "project_maps" / stem, close=True, bbox_inches="tight", pad_inches=0.02)[0]
        record_caption(out, (
            f"{_caption_row(p, in_model, dates)}. "
            + ("Footprint as the hindcast applies it (hatteras_site_config.HATTERAS_NOURISHMENT_PROJECTS). "
               if in_model else "Footprint from beach_nourishment.RECORD_ONLY_PROJECTS; dashed because the fill "
               "falls after the CoastSat data ends (2026-01-13). ")
            + f"The shaded outline is the footprint, thin white boxes the 500 m model domains, with "
            f"{MAP_PAD_DOMAINS} unfilled domains either side; the footprint's first and last domains are numbered. "
            "The volume is spread evenly over the footprint. Red lines are the Buxton groins, white lines the piers. "
            f"Imagery: {IMAGERY_NOTE}. "
            "See nourishment_summary_table for reported extents and sources."))
        outs.append(out)
    return outs


# The four model fills side by side for a paper, with a locator panel at left
def fig_project_maps_paper(model: pd.DataFrame) -> Path:
    L = _layers()
    dates = placement_dates()
    # North to south, as the locator reads top to bottom; the two Buxton fills in year order
    fills = [(p, True) for p in sorted(HATTERAS_NOURISHMENT_PROJECTS, key=lambda q: (-min(q.gis_domains), q.year))]
    # Fills on the same footprint share one panel ("Buxton, 2017 & 2022")
    panels = []
    for p, m in fills:
        same = [g for g in panels if g[0][0].gis_domains == p.gis_domains]
        (same[0] if same else panels.append([]) or panels[-1]).append((p, m))
    fs = 8.0
    # The locator (narrower, same height) then one panel per footprint; one north arrow, in the locator
    map_w, gap, left, loc_w = 0.21, 0.042, 0.02, 0.10
    box_in = map_w * PAPER_WIDTH_IN * PAPER_ASPECT              # drawn height of an equal-aspect map panel, inches
    below_in, above_in = 0.34, 0.24                             # longitude labels; titles (legend sits under the locator)
    fig_h = box_in + below_in + above_in
    fig = plt.figure(figsize=(PAPER_WIDTH_IN, fig_h))
    bottom, avail = below_in / fig_h, box_in / fig_h
    box_h, y0 = avail, bottom
    windows = []
    # Same map scale in every panel, so the latitude ticks and scale bars match
    min_h = max(L["dom"].loc[min(g[0][0].gis_domains):max(g[0][0].gis_domains)].total_bounds[3]
                - L["dom"].loc[min(g[0][0].gis_domains):max(g[0][0].gis_domains)].total_bounds[1]
                for g in panels) + 2 * PAPER_PAD_DOMAINS * SPACING_M
    for i, group in enumerate(panels):
        p, in_model = group[0]
        ax = fig.add_axes([left + loc_w + gap + i * (map_w + gap), bottom, map_w, avail])
        bx = _draw_fill(ax, p, in_model, L, PAPER_PAD_DOMAINS, PAPER_ASPECT, i, fs=fs,
                        corner=(0.03, 0.925, "left", "top"), oneline=True, north=False, title_outside=True,
                        scale_at=(0.035, 0.012) if _town(p) == "Rodanthe" else (0.825, 0.012),
                        years=" & ".join(str(q.year) for q, _ in group), min_h=min_h)
        # The locator boxes the footprint itself, not the panel window: adjacent windows (Avon, Buxton) overlap
        windows.append((f"{chr(97 + i)}", L["dom"].loc[min(p.gis_domains):max(p.gis_domains)].total_bounds))
    # The locator is shorter than the maps, top flush with them; the legend fills the space below it
    loc_frac = 0.72
    loc = fig.add_axes([left, y0 + box_h * (1 - loc_frac), loc_w, box_h * loc_frac])
    _locator(loc, L, windows, fs - 1.5, fill=True, villages=True)
    # The island named on the map itself (no panel title), lat/long on the left and top edges
    loc.apply_aspect()
    (lx0, lx1), (ly0, ly1) = loc.get_xlim(), loc.get_ylim()
    loc.text(lx0 + 0.42 * (lx1 - lx0), ly0 + 0.56 * (ly1 - ly0), "Hatteras Island", rotation=74, ha="center",
             va="center", fontsize=fs - 1.0, fontstyle="italic", color=INK, zorder=5)
    _latlon_ticks(loc, (lx0, ly0, lx1, ly1), fs - 1.5, lat_step=0.1, lon_step=0.1)
    loc.xaxis.tick_top()
    loc.tick_params(axis="x", labeltop=True, labelbottom=False)
    # A thin line north arrow, lower right of the locator over open water
    loc.annotate("", xy=(0.86, 0.13), xytext=(0.86, 0.04), xycoords="axes fraction",
                 arrowprops=dict(arrowstyle="-|>,head_length=0.5,head_width=0.18", color=INK, lw=0.7,
                                 shrinkA=0, shrinkB=0), zorder=6)
    loc.text(0.86, 0.14, "N", transform=loc.transAxes, ha="center", va="bottom", fontsize=fs - 1.5,
             color=INK, zorder=6)
    hs = _legend_handles(sorted({_town(p) for p, _ in fills}, key=["Rodanthe", "Avon", "Buxton"].index))
    hs[-1].set_label("Model domain")              # its 500 m width is in the caption; the full label overruns panel (a)
    fig.legend(handles=hs,
               loc="lower left", ncol=1, frameon=False, fontsize=fs - 1.5,
               bbox_to_anchor=(left - 0.005, y0 - 0.01), handlelength=1.4, labelspacing=0.55, borderaxespad=0.0)
    out = save(fig, OUT_DIR / "project_maps" / "nourishment_maps_paper", close=True,
               bbox_inches="tight", pad_inches=0.02)[0]
    rows = "; ".join(f"({chr(97 + i)}) " + "; ".join(_caption_row(p, m, dates) for p, m in group)
                     for i, group in enumerate(panels))
    record_caption(out, (
        "Beach nourishment projects applied in the hindcast, north up. The shaded outline in each panel is the "
        "fill footprint (colour marks the community) and thin white lines its 500 m model domains; the first and "
        "last domains are numbered. The 2017 and 2022 Buxton fills share one footprint and one panel. Each "
        "reported volume is spread evenly over its footprint. Red lines are the Buxton groins, white lines the "
        f"piers. All three panels share one scale. The left panel boxes each footprint, (a)–({chr(96 + len(panels))}), on "
        "Hatteras Island. "
        f"Imagery: {IMAGERY_NOTE}. {rows}. Reported extents and sources: nourishment_summary_table."))
    return out


# The fills as one table: CSV plus a booktabs-style figure 180 mm wide
def fig_summary_table() -> list[Path]:
    dates = placement_dates()
    rows = []
    for p, in_model in all_fills():
        words, src, url = REPORTED[(p.name, p.year)]
        s, e = dates.get(p.name, ("", ""))
        rows.append(dict(
            community=_town(p), year=p.year, placement_start=s, placement_end=e, reported_extent=words,
            source=src, source_url=url, model_gis=f"{min(p.gis_domains)}–{max(p.gis_domains)}",
            model_length_km=len(p.gis_domains) * SPACING_M / 1000, volume_yd3=int(p.volume_cubic_yards),
            volume_m3=round(p.volume_m3_total), volume_m3_per_m=round(p.volume_m3_per_m(SPACING_M), 1),
            in_model="yes" if in_model else "no"))
    tab = pd.DataFrame(rows)
    csv = OUT_DIR / "nourishment_summary_table.csv"
    tab.to_csv(csv, index=False)

    import textwrap
    cols = [("Community", 0.0, "left"), ("Placed", 0.085, "left"), ("Reported extent (source)", 0.215, "left"),
            ("Model\nGIS", 0.575, "center"), ("Length\n(km)", 0.635, "center"),
            ("Volume\n(10$^{6}$ yd$^{3}$)", 0.705, "center"), ("Volume\n(10$^{6}$ m$^{3}$)", 0.79, "center"),
            ("Density\n(" + PER_M + ")", 0.875, "center"), ("In\nmodel", 0.955, "center")]
    fs = 7.0
    body = []
    for r in tab.itertuples():
        ext = textwrap.fill(f"{r.reported_extent} ({r.source})", 58)
        placed = f"{r.placement_start}\nto {r.placement_end}"
        body.append([f"{r.community}\n{r.year}", placed, ext, r.model_gis, f"{r.model_length_km:.1f}",
                     f"{r.volume_yd3 / 1e6:.2f}", f"{r.volume_m3 / 1e6:.2f}", f"{r.volume_m3_per_m:.0f}",
                     r.in_model])
    n_lines = [max(c.count("\n") + 1 for c in row) for row in body]
    line_h = 0.135
    head_h = 2 * line_h + 0.08
    h = head_h + sum(n * line_h + 0.07 for n in n_lines) + 0.12
    fig = plt.figure(figsize=(PAPER_WIDTH_IN, h))
    ax = fig.add_axes([0.01, 0, 0.98, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(h, 0)
    ax.axis("off")
    y = 0.04
    ax.axhline(y, color=INK, lw=0.9)
    for name, x, ha in cols:
        xx = x + (0.03 if ha == "center" else 0)
        ax.text(xx, y + 0.03, name, ha=ha, va="top", fontsize=fs, fontweight="bold", color=INK, linespacing=1.15)
    y += head_h
    ax.axhline(y, color=INK, lw=0.5)
    for row, n in zip(body, n_lines):
        y += 0.035
        for (name, x, ha), cell in zip(cols, row):
            xx = x + (0.03 if ha == "center" else 0)
            ax.text(xx, y, cell, ha=ha, va="top", fontsize=fs, color=INK, linespacing=1.2)
        y += n * line_h + 0.035
    ax.axhline(y, color=INK, lw=0.9)
    ax.text(0, y + 0.03, "Density: model volume spread evenly over the footprint's 500 m domains. "
            "2026 fills fall after the CoastSat record (to 2026-01-13) and are not modelled; "
            "Buxton 2026 volume is planned, not as-built.", ha="left", va="top", fontsize=fs - 1, color=INK_MUTED)
    out = save(fig, OUT_DIR / "nourishment_summary_table", close=True, bbox_inches="tight", pad_inches=0.03)
    record_caption(out[0], (
        "Beach nourishment on the modelled reach: placement dates, the extent each source reports, and the footprint "
        "and volume the hindcast applies. Reported wording in full, with coordinates for every limit: "
        "reported_extent/reported_limits.csv; placement-date sources: datasets/nourishment_placement_dates.csv."))
    return [out[0], csv]


# Run: the figures
def main():
    apply_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    model = model_projects()
    record = record_projects()

    table = pd.concat([model, record], ignore_index=True)
    table.insert(0, "written", dt.datetime.now().strftime("%Y-%m-%d"))
    csv = OUT_DIR / "nourishment_projects.csv"
    table.to_csv(csv, index=False, float_format="%.1f")

    outs = [fig_when_where(model, record),
            fig_volume_alongshore(model),
            fig_domain_map(model)] + fig_project_maps(model) + [fig_project_maps_paper(model)] + fig_summary_table()
    print("model input:")
    print(model[["name", "year", "first_gis", "last_gis", "volume_cy",
                 "volume_m3_per_m"]].to_string(index=False))
    print("record rows:", len(record),
          "| on the reach:", int(record.in_modelled_reach.sum()),
          "| off the reach:", int((~record.in_modelled_reach).sum()))
    for o in outs + [csv]:
        print("wrote", o.relative_to(PROJECT_ROOT))


if __name__ == "__main__":
    main()
