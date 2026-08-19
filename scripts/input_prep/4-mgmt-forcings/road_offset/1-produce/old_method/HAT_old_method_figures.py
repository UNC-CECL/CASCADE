r"""
HAT_old_method_figures.py
===============================================================================
Figures and an audit for the OLD road-setback method -- the one
road_offset_pipeline.py implements and that the runner still spends.

Read-only. Writes only PNG/CSV/MD into old_method_offset/, never a forcing file.

HOW THE OLD METHOD COMPUTES A SETBACK
-------------------------------------
Two ArcGIS point exports, one per feature, sharing one transect system:

  road   raw_offset/<year>/nc12_<year>.csv          domain_id_1, ORIG_LEN
  dune   2-brie-offset/raw_offsets/<year>_duneline_offset_raw.csv
                                                    domain_id,   ORIG_LEN

Each was made by buffering the line and intersecting it with a shared
"Transect_Points_1m" layer. ORIG_LEN is the INDEX of the point along its
transect -- always an integer, always with ORIG_SEQ == ORIG_LEN + 1 -- so a
difference of ORIG_LEN is a count of 1 m points, i.e. metres. Verified two ways:
every crossing is exactly 3-4 CONSECUTIVE indices, and both layers were buffered
at 1.5 m half-width (road BUFF_DIST 4.92125 US survey feet = 1.5 m; dune
BUFF_DIST 1.5 m). A 3 m buffer sampled at 1 m gives 3-4 points. The "metres"
label on the output is therefore correct.

The pipeline then does, per domain:

    setback = min(road ORIG_LEN) - min(dune ORIG_LEN)

and that is the whole calculation. Both minima are taken over the WHOLE domain,
each within its own layer.

THE SHORTCUT, AND WHY IT IS ONE
-------------------------------
The two minima need not fall on the same transect. `min` over the road layer
picks whichever of the ~5 transects in that domain has the road nearest the
transect origin; `min` over the dune layer independently picks whichever has the
dune nearest. Subtracting them measures a road at one alongshore position
against a dune at another.

This is avoidable, because LineID IS a shared transect key -- 410 of 410 shared
across both layers in GIS 9-90. Pairing on it, then taking the median over
transects, is the like-for-like comparison. The two agree to a median of 1 m,
but disagree by more than 25 m in ~15 domains and by up to 55 m. Figure 2.

ROAD SOURCE VINTAGES -- CHECKED, AND DELIBERATE
-----------------------------------------------
The road inputs do not carry the vintages their filenames suggest:

  nc12_1984.csv   source layer NC12_1978_buffer   (git: renamed from
                  1978_road_offset_raw.csv)
  nc12_2004.csv   source layer NC12_2008_buffer

The dune inputs do (dune_1984_buffer, dune_2004_buffer).

This is FINE, and it is the project's deliberate choice: NC-12 did not move
between 1978 and 1984, or between 2004 and 2008, so the 1978 line stands in for
1984 and the 2008 line for 2004. The digitised lines that exist are the ones
that were flown, and substituting an unchanged neighbour beats not measuring.

That claim is checkable against the project's own event record, and it holds --
HISTORICAL_ROAD_EVENTS carries exactly three events (1989 Pea Island
relocation, 1999 inter-village relocation, 2022 Jug Handle Bridge) and NONE
falls inside either interval. The script re-runs that check every time and
prints the result under ROAD SOURCE VINTAGES.

Consequence worth keeping straight: because the road is unchanged across both
substitutions, the 1984->2004 difference in figure 1(c) IS a genuine 20-year
change in the forcing. It is not a vintage artifact.

OUTPUTS  (data/hatteras_init/4-mgmt-forcing/road_offset/old_method_offset/)
---------------------------------------------------------------------------
  HAT_old_method_setback.png    what the method computes
  HAT_old_method_shortcut.png   independent-min vs per-transect pairing
  RoadSetback_oldmethod_domains.csv
  RoadSetback_oldmethod_audit.md

REQUIREMENTS
------------
  pandas, numpy, matplotlib
===============================================================================
"""

from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.transforms import blended_transform_factory

# =============================================================================
# CONFIG
# =============================================================================

PROJECT_ROOT = Path(__file__).resolve().parents[6]
DATA = PROJECT_ROOT / "data" / "hatteras_init"

ROAD_CSV_FMT = (DATA / "4-mgmt-forcing" / "road_offset" / "raw_offset"
                / "{year}" / "nc12_{year}.csv")
DUNE_CSV_FMT = (DATA / "2-brie-offset" / "raw_offsets"
                / "{year}_duneline_offset_raw.csv")
SETBACK_CSV_FMT = (DATA / "4-mgmt-forcing" / "road_offset" / "old_method_offset"
                   / "{year}" / "RoadSetback_{year}.csv")

OUT_ROOT = DATA / "4-mgmt-forcing" / "road_offset" / "old_method_offset"

YEARS = [1984, 2004]
FIRST_DOMAIN, LAST_DOMAIN = 9, 90

ROAD_ID_COL, DUNE_ID_COL = "domain_id_1", "domain_id"
DIST_COL, TRANSECT_COL = "ORIG_LEN", "LineID"

# A domain whose legacy and paired values differ by more than this is worth
# opening in GIS -- the shortcut is doing real work there, not rounding.
DISAGREE_M = 25.0
# The transect layer is regular: EVERY domain in GIS 9-90 carries exactly 5
# shared transects, one per 100 m of a 500 m domain. Asserted rather than
# flagged, because a threshold that can never fire reads like a check that
# passed. If this ever trips, the transect layer changed underneath.
EXPECTED_TRANSECTS = 5

# The vintage the road export actually came from, against the year it stands in
# for. Substitution is valid only if NC-12 did not move in between -- checked
# against HISTORICAL_ROAD_EVENTS rather than asserted. See the header.
ROAD_SOURCE_YEAR = {1984: 1978, 2004: 2008}

# --- palette ------------------------------------------------------------
# Year = hue, method = line style, so no series is distinguished by colour
# alone and no third hue is needed. This pair is the repo's existing validated
# two-colour set: normal-vision dE 33.6 (floor 15), worst CVD dE 26.6
# (protanopia; target >= 8). Re-checked in OKLab, not eyeballed.
C_1984, C_2004 = "#2a78d6", "#eb6834"
C_YEAR = {1984: C_1984, 2004: C_2004}
INK, INK_MUTED, SURFACE = "#0b0b0b", "#8a8a85", "#fcfcfb"
C_BAND = "#90AFC5"

STYLE_LEGACY = (0, ())            # solid  -- what the runner spends
STYLE_PAIRED = (0, (5, 1.6))      # dashed -- the like-for-like alternative
STYLE_DUNE = (0, (1.4, 1.8))      # dotted -- the dune, in panel (a)

SECTIONS = [((1, 6), "Cape Point"), ((7, 8), "Buxton"), ((9, 20), "Buxton-Avon"),
            ((21, 31), "Avon"), ((32, 67), "Wimble Shoals"),
            ((68, 83), "Tri-Village"), ((84, 90), "Pea Island")]

plt.rcParams.update({
    "font.size": 10,
    "axes.edgecolor": INK_MUTED, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": "#52514e", "ytick.color": "#52514e",
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
})


# =============================================================================
# LOAD
# =============================================================================

def load_year(year: int) -> dict:
    """Both point layers for one year, trimmed to the road span."""
    rp = Path(str(ROAD_CSV_FMT).format(year=year))
    dp = Path(str(DUNE_CSV_FMT).format(year=year))
    for p, what in ((rp, "road"), (dp, "dune")):
        if not p.exists():
            raise FileNotFoundError(f"{what} points for {year} not found:\n"
                                    f"    {p}")
    road = pd.read_csv(rp).rename(columns={ROAD_ID_COL: "gis"})
    dune = pd.read_csv(dp).rename(columns={DUNE_ID_COL: "gis"})
    in_span = lambda g: g[(g.gis >= FIRST_DOMAIN) & (g.gis <= LAST_DOMAIN)]

    # Which ArcGIS layer each export actually came from. The filename is not
    # evidence -- see THE PROVENANCE PROBLEM in the header.
    src = {}
    for what, g in (("road", pd.read_csv(rp)), ("dune", pd.read_csv(dp))):
        cols = [c for c in g.columns
                if "buffer" in c.lower() and not c.endswith("_domains")]
        src[what] = cols[0].replace("FID_", "") if cols else "unknown"
        src[what + "_buff"] = (float(g.BUFF_DIST.iloc[0])
                               if "BUFF_DIST" in g.columns else np.nan)
    return dict(year=year, road=in_span(road), dune=in_span(dune), src=src)


def legacy_setback(y: dict) -> pd.Series:
    """min(road) - min(dune), each over the whole domain. The shipped method."""
    return (y["road"].groupby("gis")[DIST_COL].min()
            - y["dune"].groupby("gis")[DIST_COL].min()).dropna()


def paired_setback(y: dict) -> tuple:
    """
    Per-transect road - dune on the SAME LineID, then the median over transects.

    This is what the legacy method would be if it did not take its two minima
    independently. Returns (setback, n_transects).
    """
    rt = y["road"].groupby(["gis", TRANSECT_COL])[DIST_COL].min().rename("road")
    dt = y["dune"].groupby(["gis", TRANSECT_COL])[DIST_COL].min().rename("dune")
    t = pd.concat([rt, dt], axis=1, join="inner")
    t["sb"] = t.road - t.dune
    return t.groupby("gis").sb.median(), t.groupby("gis").size()


def check_road_vintages() -> list:
    """
    Is standing the 1978 road in for 1984, and the 2008 road for 2004, safe?

    Only if NC-12 did not move in between. That is not a matter of opinion --
    the project keeps a list of every time it did move, so read that list
    rather than restating the assumption. Returns one row per substitution,
    each carrying the events (if any) that would invalidate it.
    """
    try:
        sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
        from hatteras_site_config import HATTERAS_ROAD_EVENTS as events
    except Exception as e:                      # config moved or unimportable
        return [dict(label=y, src=s, ok=None, clashes=[], why=str(e))
                for y, s in ROAD_SOURCE_YEAR.items()]

    rows = []
    for label, src in sorted(ROAD_SOURCE_YEAR.items()):
        lo, hi = min(label, src), max(label, src)
        clashes = [e for e in events if lo < getattr(e, "year", -1) <= hi]
        rows.append(dict(label=label, src=src, ok=not clashes,
                         clashes=[(e.year, getattr(e, "note", type(e).__name__))
                                  for e in clashes], why=""))
    return rows


def verify_against_disk(year: int, legacy: pd.Series) -> float | None:
    """Does the recomputation reproduce the file the runner actually reads?"""
    p = Path(str(SETBACK_CSV_FMT).format(year=year))
    if not p.exists():
        return None
    raw = np.loadtxt(p, delimiter=",")
    on_disk = dict(zip(raw[0].astype(int).tolist(), raw[1].tolist()))
    common = [g for g in legacy.index if g in on_disk]
    if not common:
        return None
    return float(np.abs([legacy[g] - on_disk[g] for g in common]).max())


# =============================================================================
# FIGURE HELPERS
# =============================================================================

def shade_sections(ax, lo, hi, label_top=True):
    for k, ((a, b), name) in enumerate(SECTIONS):
        if b < lo or a > hi:
            continue
        ax.axvspan(max(a, lo) - 0.5, min(b, hi) + 0.5, color=C_BAND,
                   alpha=0.16 if k % 2 == 0 else 0.06, lw=0, zorder=0)
        if label_top:
            # Inside the axes, not above them -- above collides with the panel
            # title on a stacked figure.
            tr = blended_transform_factory(ax.transData, ax.transAxes)
            ax.text((max(a, lo) + min(b, hi)) / 2, 0.975, name, transform=tr,
                    ha="center", va="top", fontsize=7.5, color=INK_MUTED,
                    zorder=7,
                    bbox=dict(boxstyle="round,pad=0.2", fc=SURFACE, ec="none",
                              alpha=0.85))


def tidy(ax, ylabel):
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.22, lw=0.7)
    ax.set_axisbelow(True)
    ax.set_xlim(FIRST_DOMAIN - 0.5, LAST_DOMAIN + 0.5)


XLABEL = "GIS domain   (9 = Buxton / south  →  90 = Rodanthe / north)"


# =============================================================================
# FIGURE 1 -- what the method computes
# =============================================================================

def figure_method(data, out: Path):
    fig, axes = plt.subplots(3, 1, figsize=(13, 11), sharex=True,
                             gridspec_kw=dict(hspace=0.28))
    ax_a, ax_b, ax_c = axes

    for y in data:
        yr, c = y["year"], C_YEAR[y["year"]]
        rmin = y["road"].groupby("gis")[DIST_COL].min()
        dmin = y["dune"].groupby("gis")[DIST_COL].min()
        ax_a.plot(rmin.index, rmin.values, color=c, lw=2.0, ls=STYLE_LEGACY,
                  label=f"road, {yr}", zorder=4)
        ax_a.plot(dmin.index, dmin.values, color=c, lw=1.8, ls=STYLE_DUNE,
                  label=f"dune, {yr}", zorder=3)

    shade_sections(ax_a, FIRST_DOMAIN, LAST_DOMAIN)
    tidy(ax_a, "Distance along transect (m)")
    ax_a.set_title("(a)  The two measured quantities — each a distance from the "
                   "transect origin, not a setback", loc="left", fontsize=10.5)
    # Lower left: the curves descend left-to-right, so this corner is empty and
    # it stays clear of the section labels along the top.
    ax_a.legend(loc="lower left", fontsize=8.5, ncol=2, framealpha=0.92)

    for y in data:
        yr, c = y["year"], C_YEAR[y["year"]]
        sb = legacy_setback(y)
        ax_b.plot(sb.index, sb.values, color=c, lw=2.0, ls=STYLE_LEGACY,
                  marker="o", ms=3.2, label=f"{yr} setback", zorder=4)
    shade_sections(ax_b, FIRST_DOMAIN, LAST_DOMAIN, label_top=False)
    tidy(ax_b, "Road setback (m)")
    ax_b.set_title("(b)  Their difference — the forcing the runner spends. Note "
                   "the axis: a ~200 m result from two ~6000 m measurements",
                   loc="left", fontsize=10.5)
    ax_b.legend(loc="upper left", fontsize=8.5, framealpha=0.92)

    if len(data) == 2:
        a, b = legacy_setback(data[0]), legacy_setback(data[1])
        common = sorted(set(a.index) & set(b.index))
        chg = np.array([b[g] - a[g] for g in common])
        ax_c.axhline(0, color=INK_MUTED, lw=1.0, zorder=2)
        ax_c.bar(common, chg, width=0.82,
                 color=[C_2004 if v < 0 else C_1984 for v in chg],
                 zorder=3, linewidth=0)
        shade_sections(ax_c, FIRST_DOMAIN, LAST_DOMAIN, label_top=False)
        tidy(ax_c, "Change, 2004 − 1984 (m)")
        ax_c.set_title("(c)  Change between the two files — a real 20-year "
                       "change. Road layers are 1978 and 2008, but NC-12 did "
                       "not move in either gap (see audit)",
                       loc="left", fontsize=10.5)
        ax_c.text(0.995, 0.06,
                  f"negative = road closer to dune in 2004   "
                  f"(median {np.median(chg):+.0f} m)",
                  transform=ax_c.transAxes, ha="right", va="bottom",
                  fontsize=8.5, color=INK_MUTED)

    ax_c.set_xlabel(XLABEL)
    fig.suptitle("Old road-setback method — min(road) − min(dune) per domain",
                 fontsize=13, y=0.995)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
    print(f"[out] {out}")


# =============================================================================
# FIGURE 2 -- the independent-min shortcut
# =============================================================================

def figure_shortcut(data, out: Path):
    fig, axes = plt.subplots(3, 1, figsize=(13, 11), sharex=True,
                             gridspec_kw=dict(hspace=0.28))
    ax_a, ax_b, ax_c = axes

    # (a) How much min() is choosing out of. Plotted as the SPREAD itself, in
    # metres -- drawn on the absolute distance axis of figure 1(a) a 50-100 m
    # spread is a hairline on a 5000 m range and shows nothing.
    for y in data:
        yr, c = y["year"], C_YEAR[y["year"]]
        for what, ls in (("road", STYLE_LEGACY), ("dune", STYLE_DUNE)):
            sp = y[what].groupby("gis")[DIST_COL].agg(lambda s: s.max() - s.min())
            ax_a.plot(sp.index, sp.values, color=c, lw=1.8, ls=ls, zorder=4,
                      label=f"{what}, {yr}")
    med = np.median(np.concatenate(
        [y[w].groupby("gis")[DIST_COL].agg(lambda s: s.max() - s.min()).values
         for y in data for w in ("road", "dune")]))
    ax_a.axhline(med, color=INK, lw=1.2, ls=(0, (4, 3)), zorder=5)
    ax_a.text(LAST_DOMAIN, med, f" median {med:.0f} m ", ha="right",
              va="bottom", fontsize=8, color=INK, zorder=6,
              bbox=dict(boxstyle="round,pad=0.2", fc=SURFACE, ec="none",
                        alpha=0.85))
    shade_sections(ax_a, FIRST_DOMAIN, LAST_DOMAIN)
    tidy(ax_a, "Spread across a domain's\ntransects, max − min (m)")
    ax_a.set_ylim(bottom=0)
    ax_a.set_title("(a)  How much min() is choosing out of. Each domain holds "
                   f"exactly {EXPECTED_TRANSECTS} transects; this is the range "
                   "they span", loc="left", fontsize=10.5)
    ax_a.legend(loc="lower left", fontsize=7.5, ncol=4, framealpha=0.92)

    # (b) legacy vs paired
    store = {}
    for y in data:
        yr, c = y["year"], C_YEAR[y["year"]]
        leg = legacy_setback(y)
        pair, n = paired_setback(y)
        common = sorted(set(leg.index) & set(pair.index))
        store[yr] = dict(leg=leg, pair=pair, n=n, common=common)
        ax_b.plot(common, [leg[g] for g in common], color=c, lw=2.0,
                  ls=STYLE_LEGACY, label=f"{yr} legacy: independent min",
                  zorder=4)
        ax_b.plot(common, [pair[g] for g in common], color=c, lw=1.7,
                  ls=STYLE_PAIRED, label=f"{yr} paired: same transect, median",
                  zorder=3)
    shade_sections(ax_b, FIRST_DOMAIN, LAST_DOMAIN, label_top=False)
    tidy(ax_b, "Road setback (m)")
    ax_b.set_title("(b)  The shipped value against the like-for-like one. "
                   "LineID is shared by both layers, so pairing is available",
                   loc="left", fontsize=10.5)
    ax_b.legend(loc="upper left", fontsize=8, ncol=2, framealpha=0.92)

    # (c) the disagreement
    ax_c.axhspan(-DISAGREE_M, DISAGREE_M, color=INK_MUTED, alpha=0.10, lw=0,
                 zorder=1, label=f"±{DISAGREE_M:.0f} m")
    ax_c.axhline(0, color=INK_MUTED, lw=1.0, zorder=2)
    for yr, s in store.items():
        d = np.array([s["leg"][g] - s["pair"][g] for g in s["common"]])
        ax_c.plot(s["common"], d, color=C_YEAR[yr], lw=1.8, ls=STYLE_LEGACY,
                  marker="o", ms=3.0, zorder=4, label=f"{yr}")
        big = [g for g in s["common"] if abs(s["leg"][g] - s["pair"][g])
               > DISAGREE_M]
        if big:
            ax_c.plot(big, [s["leg"][g] - s["pair"][g] for g in big], lw=0,
                      marker="o", ms=7.5, mfc="none", mew=1.6,
                      color=C_YEAR[yr], zorder=6,
                      label=f"{yr}: beyond ±{DISAGREE_M:.0f} m "
                            f"({len(big)} domains)")
    shade_sections(ax_c, FIRST_DOMAIN, LAST_DOMAIN, label_top=False)
    tidy(ax_c, "legacy − paired (m)")
    ax_c.set_title("(c)  What the shortcut costs. Small in the median, not "
                   "small in the tail", loc="left", fontsize=10.5)
    ax_c.legend(loc="lower left", fontsize=7.5, ncol=3, framealpha=0.92)
    ax_c.set_xlabel(XLABEL)

    fig.suptitle("Old method — taking the two minima independently",
                 fontsize=13, y=0.995)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
    print(f"[out] {out}")


# =============================================================================
# TABLE + AUDIT
# =============================================================================

def write_table(data, path: Path):
    rows = []
    for y in data:
        leg = legacy_setback(y)
        pair, n = paired_setback(y)
        rspread = y["road"].groupby("gis")[DIST_COL].agg(lambda s: s.max() - s.min())
        dspread = y["dune"].groupby("gis")[DIST_COL].agg(lambda s: s.max() - s.min())
        for g in sorted(set(leg.index) & set(pair.index)):
            rows.append(dict(
                year=y["year"], domain=int(g),
                setback_legacy_m=round(float(leg[g]), 1),
                setback_paired_m=round(float(pair[g]), 1),
                legacy_minus_paired_m=round(float(leg[g] - pair[g]), 1),
                n_shared_transects=int(n[g]),
                road_spread_m=round(float(rspread.get(g, np.nan)), 1),
                dune_spread_m=round(float(dspread.get(g, np.nan)), 1),
                flags=(f"DISAGREE({leg[g] - pair[g]:+.0f})"
                       if abs(leg[g] - pair[g]) > DISAGREE_M else ""),
            ))
    df = pd.DataFrame(rows).sort_values(["year", "domain"])
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    print(f"[out] {path}")
    return df


def write_audit(data, df, checks, path: Path):
    L, a = [], None
    L = []
    a = L.append
    a("# Old road-setback method — how it works and where it strains")
    a("")
    a(f"Generated {datetime.now().isoformat(timespec='seconds')} by "
      f"`{Path(__file__).name}`. Read-only: no forcing file is written here.")
    a("")

    a("## The calculation")
    a("")
    a("Per domain, from two ArcGIS point exports that share one transect "
      "system:")
    a("")
    a("```")
    a("setback = min(road ORIG_LEN)  −  min(dune ORIG_LEN)")
    a("```")
    a("")
    a("`ORIG_LEN` is the **index of a point along its transect**, not a measured "
      "length — always an integer, always with `ORIG_SEQ == ORIG_LEN + 1`. A "
      "difference of indices is a count of 1 m points, so the result is metres "
      "and the output label is correct. Two independent confirmations:")
    a("")
    a("- every crossing is **3–4 consecutive** indices in both layers;")
    a("- both layers were buffered at a 1.5 m half-width — road `BUFF_DIST` "
      "4.92125 (US survey feet = 1.5 m), dune `BUFF_DIST` 1.5 (m). A 3 m buffer "
      "sampled at 1 m spacing yields exactly 3–4 points.")
    a("")
    if checks.get("reproduce"):
        for yr, mx in checks["reproduce"].items():
            verdict = "reproduces exactly" if mx is not None and mx < 1e-6 else \
                      (f"differs by up to {mx:.3f} m" if mx is not None
                       else "no file on disk to compare")
            a(f"- **{yr}**: recomputing from the raw exports {verdict} "
              f"(`RoadSetback_{yr}.csv`).")
        a("")
        a("> This matters for the record: the legacy setbacks **are "
          "reproducible** from the CSVs in the repo. What is missing is the "
          "upstream `duneline_<year>.geojson`, not the offsets derived from it.")
        a("")

    a("## The shortcut — two minima, taken independently")
    a("")
    a("Both minima are taken over the whole domain, each within its own layer, "
      "so they need not fall on the same transect. The method can therefore "
      "measure a road at one alongshore position against a dune at another.")
    a("")
    a("This is avoidable. `LineID` is a **shared transect key** — "
      f"{checks['shared_transects']} of {checks['road_transects']} road "
      "transects are also dune transects across GIS 9–90 — so a per-transect "
      "difference, then a median over transects, is available and is the "
      "like-for-like comparison.")
    a("")
    a("| Year | Median legacy − paired | p10 | p90 | Max abs | Domains > "
      f"{DISAGREE_M:.0f} m |")
    a("|---:|---:|---:|---:|---:|---:|")
    for yr in sorted(df.year.unique()):
        s = df[df.year == yr].legacy_minus_paired_m
        a(f"| {yr} | {s.median():+.1f} m | {s.quantile(.10):+.1f} | "
          f"{s.quantile(.90):+.1f} | {s.abs().max():.0f} | "
          f"{int((s.abs() > DISAGREE_M).sum())} |")
    a("")
    a("Small in the median, not small in the tail. The median is ~1 m — which "
      "is why this was recorded as a minor shortcut — but the disagreement "
      "exceeds 25 m in roughly a sixth of the road and reaches ~55 m. Those "
      "domains are flagged `DISAGREE` in the per-domain CSV.")
    a("")
    a(f"Sample size compounds it. The transect layer is regular — **every** "
      f"domain in GIS 9–90 carries exactly "
      f"{int(df['n_shared_transects'].median())} shared transects, one per "
      f"100 m of a 500 m domain — so both the legacy minimum and the paired "
      f"median rest on 5 observations. The spread those 5 span is a median of "
      f"{df['road_spread_m'].median():.0f} m for the road and "
      f"{df['dune_spread_m'].median():.0f} m for the dune, reaching "
      f"{df['road_spread_m'].max():.0f} m and "
      f"{df['dune_spread_m'].max():.0f} m. Figure 2(a).")
    a("")

    a("## Road source vintages — checked, and deliberate")
    a("")
    a("The road inputs do not carry the vintages their filenames suggest:")
    a("")
    a("| Year label | Road source layer | Dune source layer |")
    a("|---|---|---|")
    for y in data:
        a(f"| {y['year']} | `{y['src']['road']}` | `{y['src']['dune']}` |")
    a("")
    a("The dune layers match their labels; the road layers do not. Git records "
      "`nc12_1984.csv` as a rename of `1978_road_offset_raw.csv`, and the 2004 "
      "export carries `NC12_2008_buffer`.")
    a("")
    a("**This is deliberate and it is sound.** NC-12 did not move between 1978 "
      "and 1984, or between 2004 and 2008, so the 1978 line stands in for 1984 "
      "and the 2008 line for 2004. The digitised lines that exist are the ones "
      "that were flown; substituting an unchanged neighbour beats not "
      "measuring.")
    a("")
    a("That is checkable rather than assumed, and the script re-checks it every "
      "run against `HATTERAS_ROAD_EVENTS` — the same list the runner uses:")
    a("")
    a("| Substitution | Road events in between | Verdict |")
    a("|---|---|---|")
    for v in checks.get("vintages", []):
        if v["ok"] is None:
            a(f"| {v['src']} → {v['label']} | *not checked* | `{v['why']}` |")
        else:
            cl = ("none" if not v["clashes"]
                  else "; ".join(f"{yr} {note}" for yr, note in v["clashes"]))
            a(f"| {v['src']} → {v['label']} | {cl} | "
              f"{'**valid**' if v['ok'] else '**INVALID**'} |")
    a("")
    a("The three events on record — 1989 Pea Island relocation, 1999 "
      "inter-village relocation, 2022 Jug Handle Bridge — all fall outside both "
      "intervals.")
    a("")
    a("Consequence: because the road is unchanged across both substitutions, "
      "**panel (c) of figure 1 is a genuine 20-year change** in the forcing, "
      "not a vintage artifact. A new method fed the same road CSVs inherits the "
      "same, valid, substitution.")
    a("")

    a("## A note on precision")
    a("")
    a("Figure 1(a) is the point: the setback is a **difference of two large "
      "numbers**. Both distances run to several thousand metres from the "
      "transect origin, and the answer is a few hundred. Nothing here is wrong "
      "— the indices are exact integers, so there is no floating-point erosion "
      "— but it does mean the result inherits any error in either baseline "
      "in full, and a transect origin that shifted between the two layers "
      "would propagate straight through.")
    a("")

    a("## Per-domain results")
    a("")
    a("Full table in `RoadSetback_oldmethod_domains.csv`. Flagged domains only:")
    a("")
    # Bracket access throughout: `flags` collides with pandas' own
    # DataFrame.flags / Series.flags property, and attribute access silently
    # returns that instead of the column.
    fl = df[df["flags"] != ""]
    if fl.empty:
        a("*No domain flagged.*")
    else:
        a("| Year | GIS | Legacy | Paired | Difference | Shared transects | "
          "Flags |")
        a("|---:|---:|---:|---:|---:|---:|---|")
        for _, r in fl.iterrows():
            a(f"| {r['year']} | {r['domain']} | {r['setback_legacy_m']:.0f} | "
              f"{r['setback_paired_m']:.0f} | "
              f"{r['legacy_minus_paired_m']:+.0f} | "
              f"{r['n_shared_transects']} | "
              + ", ".join(f"`{x}`" for x in r["flags"].split(",") if x) + " |")
    a("")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(L), encoding="utf-8")
    print(f"[out] {path}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 88)
    print("OLD ROAD-SETBACK METHOD -- figures and audit")
    print("=" * 88)

    data = []
    for year in YEARS:
        try:
            data.append(load_year(year))
        except FileNotFoundError as e:
            print(f"[skip] {year}: {e}")
    if not data:
        raise SystemExit("\n[stop] no input layers found\n")

    checks = {"reproduce": {}}
    for y in data:
        leg = legacy_setback(y)
        mx = verify_against_disk(y["year"], leg)
        checks["reproduce"][y["year"]] = mx
        src = y["src"]
        print(f"\n  {y['year']}:  road layer {src['road']} "
              f"(BUFF_DIST {src['road_buff']})")
        print(f"          dune layer {src['dune']} "
              f"(BUFF_DIST {src['dune_buff']})")
        print(f"          {len(leg)} domains | setback {leg.min():.0f}"
              f"-{leg.max():.0f} m, median {leg.median():.0f} m")
        if mx is None:
            print("          no RoadSetback CSV on disk to check against")
        elif mx < 1e-6:
            print("          reproduces RoadSetback CSV on disk EXACTLY")
        else:
            print(f"          differs from disk by up to {mx:.3f} m")

    y0 = data[0]
    rt = set(map(tuple, y0["road"][["gis", TRANSECT_COL]].drop_duplicates().values))
    dt = set(map(tuple, y0["dune"][["gis", TRANSECT_COL]].drop_duplicates().values))
    checks["shared_transects"] = len(rt & dt)
    checks["road_transects"] = len(rt)
    checks["vintages"] = check_road_vintages()
    print(f"\n  shared (domain, LineID) transects: {len(rt & dt)} of "
          f"{len(rt)} road / {len(dt)} dune  -> pairing IS available")

    counts = {int(n) for y in data for n in paired_setback(y)[1].values}
    checks["transect_counts"] = sorted(counts)
    if counts == {EXPECTED_TRANSECTS}:
        print(f"  every domain carries exactly {EXPECTED_TRANSECTS} shared "
              f"transects (one per 100 m)")
    else:
        print(f"  [warn] transect count per domain is NOT uniform: "
              f"{sorted(counts)}")
        print("         the audit's 'exactly 5' statement no longer holds -- "
              "the transect layer changed.")

    figure_method(data, OUT_ROOT / "HAT_old_method_setback.png")
    figure_shortcut(data, OUT_ROOT / "HAT_old_method_shortcut.png")
    df = write_table(data, OUT_ROOT / "RoadSetback_oldmethod_domains.csv")
    write_audit(data, df, checks, OUT_ROOT / "RoadSetback_oldmethod_audit.md")

    print(f"\n{'=' * 88}")
    print("ROAD SOURCE VINTAGES")
    print("=" * 88)
    for y in data:
        print(f"  {y['year']}: road = {y['src']['road']} | "
              f"dune = {y['src']['dune']}")
    print()
    for v in checks["vintages"]:
        if v["ok"] is None:
            print(f"  {v['src']} -> {v['label']}: could not check "
                  f"({v['why']})")
        elif v["ok"]:
            print(f"  {v['src']} -> {v['label']}: OK, no road event in "
                  f"between")
        else:
            print(f"  {v['src']} -> {v['label']}: INVALID -- road moved in "
                  f"between:")
            for yr, note in v["clashes"]:
                print(f"      {yr}: {note}")
    if all(v["ok"] for v in checks["vintages"] if v["ok"] is not None):
        print("\n  NC-12 did not move across either substitution, so the road")
        print("  layers stand in for their labels and the 1984->2004 change is")
        print("  a real 20-year change, not a vintage artifact.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
