"""
beach_nourishment.py
==============================================================================
When and where the beach was nourished: figures for
data/hatteras_init/4-mgmt-forcing/nourishment/.

THE QUESTION
    The hindcast fires the projects in HATTERAS_NOURISHMENT_PROJECTS
    (scripts/site_layer/hatteras_site_config.py) when they fall inside a run window.
    That list is three entries and the record behind it is a spreadsheet
    (Hatteras_Management_Timelines.xlsx, sheet Nourishment_Timeline), and
    nothing in the data tree shows the two side by side, or shows the reader
    which stretch of the island was filled in which year. These figures do.

TWO SOURCES, DRAWN TOGETHER, NEVER MERGED
    * the MODEL INPUT: the site-config projects, their extents and their
      volumes spread to m^3/m over 500 m domains. This is what a run receives.
    * the RECORD: the spreadsheet. Its domain flags for the same three
      projects are narrower than the site config's (Rodanthe 85-88 against
      84-89; Avon 23-26 against 21-28; Buxton the same 6-15) because the site
      config re-derived the footprints from the project descriptions -- the
      reasons are in the comments beside each entry. The record also carries
      fourteen Pea Island / Oregon Inlet navigation fills (1990-2004, 2013)
      that lie NORTH of GIS 90, off the modelled reach, and that no run sees.
    The figures draw the model extent as the fill and the record's flags as
    a darker inner bar, so a difference is visible rather than reconciled.

OUTPUTS   data/hatteras_init/4-mgmt-forcing/nourishment/
    nourishment_when_where.png/.pdf       year x domain event chart
    nourishment_volume_alongshore.png/.pdf m^3/m per domain, as delivered
    nourishment_domain_map.png/.pdf       the filled domains on the island
    nourishment_projects.csv              every row drawn, both sources
    CAPTIONS.md                           the text that is not on the canvas

RUN
    python scripts/input_prep/4-mgmt-forcings/beach_nourishment.py
==============================================================================
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
    save, town_bands,
)
from site_layer.hatteras_site_config import (  # noqa: E402
    HATTERAS_ANNOTATIONS, HATTERAS_DOMAINS, HATTERAS_NOURISHMENT_PROJECTS,
    HATTERAS_PERIODS,
)

DATA_DIR = PROJECT_ROOT / "data" / "hatteras_init"
from site_layer.hat_topo_version import MGMT_ROOT as MGMT_DIR, MGMT_RECORD_XLSX as RECORD_XLSX  # noqa: E402
from site_layer.hat_observed_rates import DOMAIN_BOXES as DOMAIN_FILE  # noqa: E402
from site_layer.hat_map_layers import ISLAND_OUTLINE as OUTLINE_FILE  # noqa: E402
OUT_DIR = MGMT_DIR / "nourishment"

SPACING_M = HATTERAS_DOMAINS.domain_spacing_m
N_DOMAINS = HATTERAS_DOMAINS.num_real_domains
YEAR_LO, YEAR_HI = 1984, 2024

# Vintage colours: the earlier fill is red, the later one blue, as everywhere
# two vintages share a figure (hat_figure_style).
YEAR_COLOUR = {2014: (C_1984, C_1984_FILL), 2022: (C_1997, C_1997_FILL)}


# =============================================================================
# THE TWO SOURCES
# =============================================================================

def model_projects() -> pd.DataFrame:
    """The site-config list, one row per project, with the volumes a run gets."""
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


def record_projects() -> pd.DataFrame:
    """The spreadsheet: one row per year that carries a note or a flag.
    Domain flags (1) give the footprint; a note with no flags is a fill the
    record places outside the 90 domains (Pea Island / Oregon Inlet)."""
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
        # One year can flag two separate stretches (2022: Buxton AND Avon),
        # so each contiguous run of flags is its own row.
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


def _contiguous_runs(numbers):
    runs = []
    for n in sorted(numbers):
        if runs and n == runs[-1][-1] + 1:
            runs[-1].append(n)
        else:
            runs.append([n])
    return runs


def period_windows():
    """(start, end) of every hindcast window in the site config."""
    return sorted((int(k), int(v["end_year"])) for k, v in HATTERAS_PERIODS.items())


# =============================================================================
# FIGURE 1: WHEN x WHERE
# =============================================================================

def fig_when_where(model: pd.DataFrame, record: pd.DataFrame) -> Path:
    fig, ax = plt.subplots(figsize=figsize("double", aspect=0.52),
                           constrained_layout=True)
    windows = period_windows()
    # Left margin: the hindcast windows as vertical brackets. Right margin: the
    # off-reach fills. Both are OFF the domain axis, which runs 1-90.
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

    handles = [
        Patch(facecolor=C_1984_FILL, edgecolor=C_1984, lw=0.7, label="2014 fill, model input extent"),
        Patch(facecolor=C_1997_FILL, edgecolor=C_1997, lw=0.7, label="2022 fills, model input extent"),
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
    caption(fig, (
        f"Beach nourishment on the modelled reach, {YEAR_LO}–{YEAR_HI}, by year placed (vertical) "
        f"and GIS domain (horizontal, south to north, {SPACING_M:.0f} m each). Each filled bar is one "
        f"project as the hindcast receives it from hatteras_site_config.HATTERAS_NOURISHMENT_PROJECTS: "
        + "; ".join(f"{n} {int(row.year)}, domains {int(row.first_gis)}–{int(row.last_gis)} "
                    f"({row.n_domains} domains, {row.length_km:.1f} km, {row.volume_cy/1e6:.1f} M cy)"
                    for n, row in m.iterrows())
        + ". Red is the earlier fill, blue the later ones. The dark inner bar is the footprint flagged for "
          "the same project in Hatteras_Management_Timelines.xlsx (Nourishment_Timeline sheet): "
        + "; ".join(f"{int(y)} domains {int(row.first_gis)}–{int(row.last_gis)}" for y, row in r.iterrows())
        + ". Where the two differ the site-config extent was re-derived from the project description "
          "(Rodanthe: 2 mi north of the village, stopping short of the locked boundary domain 90; "
          "Avon: Due East Road to Askins Creek North Drive) and the record's flags are kept as drawn, "
          "not reconciled. Hollow markers at the right are the {n_off} Pea Island / Oregon Inlet "
          f"navigation and emergency fills the record lists ({off_years}); they lie north of domain 90 "
          f"and no run sees them. Brackets at the left are the {len(windows)} hindcast windows; a run "
          "fires whatever falls inside its window, so "
        + _windows_sentence(windows, model.year.tolist())
        + ". Named bands along the foot are the villages."
    ).replace("{n_off}", str(n_off)))
    return save(fig, OUT_DIR / "nourishment_when_where", close=True)[0]


def _windows_sentence(windows, fill_years):
    """'1984–2004 and 1996–2010 carry no fill, 2004–2024 and 2010–2024 carry
    all three', computed from the config so the caption cannot go stale."""
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


# =============================================================================
# FIGURE 2: VOLUME ALONGSHORE
# =============================================================================

def fig_volume_alongshore(model: pd.DataFrame) -> Path:
    fig, ax = plt.subplots(figsize=figsize("double", aspect=0.40),
                           constrained_layout=True)
    per_domain = np.zeros(N_DOMAINS + 1)
    for p in HATTERAS_NOURISHMENT_PROJECTS:
        v = p.volume_m3_per_m(SPACING_M)
        dark, light = YEAR_COLOUR[p.year]
        gis = np.array(p.gis_domains)
        ax.bar(gis, [v] * len(gis), width=0.86, facecolor=light, edgecolor=dark,
               lw=0.6, zorder=3)
        per_domain[gis] += v
        ax.text(gis.mean(), v + 12, f"{p.name.split()[0]} {p.year}\n{v:.0f} m³/m",
                ha="center", va="bottom", fontsize=7.5, color=INK, zorder=6)

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
        ax.text(x + 0.4, per_domain.max() * 0.92, name, rotation=90, ha="left",
                va="top", fontsize=6.5, color=INK_MUTED)
    for name, (x, _) in HATTERAS_ANNOTATIONS.piers.items():
        ax.plot([x], [0], marker="^", ms=5, color=INK, clip_on=False, zorder=7)
        ax.text(x + 0.4, per_domain.max() * 0.04, name, rotation=90, ha="left",
                va="bottom", fontsize=6.5, color=INK_MUTED, zorder=7)

    handles = [
        Patch(facecolor=C_1984_FILL, edgecolor=C_1984, lw=0.6, label="placed 2014"),
        Patch(facecolor=C_1997_FILL, edgecolor=C_1997, lw=0.6, label="placed 2022"),
        Line2D([], [], color=INK, lw=0.8, ls=(0, (2, 2)), label="groin"),
        Line2D([], [], marker="^", ls="none", color=INK, markersize=5, label="pier"),
    ]
    fig.legend(handles=handles, loc="outside lower center", ncol=4, frameon=False)

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
        + ". Even spreading is an assumption; real fill templates taper at their ends. No domain is "
          "filled twice, so the bars are also the 1984–2024 cumulative. The 2014 Rodanthe fill is "
          "about twice as dense as the 2022 fills because it put more sand into a shorter reach. "
          "The dotted line is the Buxton groin field, from which the Buxton footprint was measured "
          "north; triangles are the piers, and the Avon footprint was placed about the pier. Named "
          "bands along the top are the villages: the Buxton and Rodanthe footprints extend out of "
          "their villages into the road corridor, the Avon footprint stays inside its village."
    ))
    return save(fig, OUT_DIR / "nourishment_volume_alongshore", close=True)[0]


# =============================================================================
# FIGURE 3: THE DOMAIN MAP
# =============================================================================

def fig_domain_map(model: pd.DataFrame) -> Path:
    import geopandas as gpd
    from shapely.affinity import rotate as shapely_rotate
    from shapely.ops import unary_union

    domains = gpd.read_file(DOMAIN_FILE)
    outline = gpd.read_file(OUTLINE_FILE).to_crs(domains.crs)
    origin = tuple(unary_union(domains.geometry.tolist()).centroid.coords[0])

    def to_strip(frame):
        # Rotate 90 deg clockwise so south is left, north right, ocean below;
        # rotation preserves distance, so the scale bar holds.
        return frame.set_geometry(frame.geometry.apply(
            lambda g: shapely_rotate(g, -90, origin=origin)), crs=frame.crs)

    strip = to_strip(domains)
    strip_outline = to_strip(outline)
    year_of = {}
    for p in HATTERAS_NOURISHMENT_PROJECTS:
        for g in p.gis_domains:
            year_of[g] = p.year
    strip["fill_year"] = strip["domain_id"].map(year_of)

    fig, ax = plt.subplots(figsize=figsize("double", aspect=0.36),
                           constrained_layout=True)
    minx, miny, maxx, maxy = strip.total_bounds
    strip_outline.plot(ax=ax, facecolor="0.93", edgecolor="none", zorder=0)
    strip.plot(ax=ax, facecolor="none", edgecolor="0.6", linewidth=0.35, zorder=1)
    for year, (dark, light) in YEAR_COLOUR.items():
        sub = strip[strip.fill_year == year]
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
    # Every tenth domain, unless a footprint end sits within one domain of it
    # (20 beside 21, 90 beside 89), plus the footprint ends themselves.
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

    for _, p in model.iterrows():
        b = strip[strip.domain_id.between(p.first_gis, p.last_gis)].total_bounds
        ax.annotate("", xy=(b[0], maxy + 0.30 * dy), xytext=(b[2], maxy + 0.30 * dy),
                    arrowprops=dict(arrowstyle="|-|", linewidth=1.0, color=INK,
                                    mutation_scale=3), annotation_clip=False)
        ax.text((b[0] + b[2]) / 2, maxy + 0.34 * dy,
                f"{p['name'].split()[0]} {int(p.year)}\ndomains {int(p.first_gis)}–{int(p.last_gis)}",
                ha="center", va="bottom", fontsize=8, color=INK)

    ax.text(minx, miny - 0.16 * dy, "Cape Point (south)", ha="left", va="top",
            fontsize=8.5, color=INK_MUTED)
    ax.text(maxx, miny - 0.16 * dy,
            "Pea Island (north)\nOregon Inlet fills lie beyond this end, not modelled",
            ha="right", va="top", fontsize=8.5, color=INK_MUTED, linespacing=1.3)
    ax.text((minx + maxx) / 2, miny - 0.10 * dy, "Atlantic Ocean", ha="center",
            va="top", fontsize=8.5, color=INK_MUTED, style="italic")

    # Scale bar in data units, 5 km, and a north arrow that points along the
    # strip since the map is rotated off north.
    bx, by = minx + 0.02 * (maxx - minx), miny - 0.40 * dy
    ax.plot([bx, bx + 5000], [by, by], color=INK, lw=2.2, solid_capstyle="butt", zorder=12)
    ax.text(bx + 2500, by + 0.05 * dy, "5 km", ha="center", va="bottom", fontsize=8, color=INK)
    # Beside the scale bar, clear of the north-end labels.
    ax0, ax1 = bx + 7500, bx + 10000
    ax.annotate("", xy=(ax1, by), xytext=(ax0, by),
                arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.0, mutation_scale=11))
    ax.text((ax0 + ax1) / 2, by + 0.05 * dy, "N", ha="center", va="bottom",
            fontsize=8.5, fontweight="bold", color=INK)

    handles = [
        Patch(facecolor=C_1984_FILL, edgecolor=C_1984, lw=0.9, label="nourished 2014"),
        Patch(facecolor=C_1997_FILL, edgecolor=C_1997, lw=0.9, label="nourished 2022"),
        Patch(facecolor="white", edgecolor="0.6", lw=0.35, label="domain, never nourished"),
    ]
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False)

    n_filled = len(year_of)
    caption(fig, (
        f"The {len(domains)} GIS domains in place, the island rotated 90° clockwise so that south is "
        "to the left, north to the right and the Atlantic at the bottom; rotation preserves distance, "
        f"so the scale bar holds. The {n_filled} domains that receive a nourishment project in the "
        "hindcast are filled by the year it was placed, with the site-config extents of Figure 1; "
        "grey is the island outline and white domains are never nourished. Every tenth domain and "
        "the ends of each footprint are numbered along the top. The Pea Island / Oregon Inlet fills "
        "in the management record lie beyond the north end of the strip."
    ))
    return save(fig, OUT_DIR / "nourishment_domain_map", close=True)[0]


# =============================================================================
# MAIN
# =============================================================================

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
            fig_domain_map(model)]
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
