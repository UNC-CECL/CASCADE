r"""
coastsat_vs_duneline.py
==============================================================================
Does the digitized dune line move with the CoastSat shoreline?

WHAT IS COMPARED  (NET CHANGE ON BOTH SIDES since 2026-09-18; Hannah: "I
wanted the coastsat vs duneline comparison to both be using net position
change")
    Dune side      3-rates/duneline/endpoint/<window>/: the end dune line minus
                   the start line per 100 m transect, domain mean, seaward
                   positive.
    Shoreline side 3-rates/coastsat/endpoint/<window>/: the mean CoastSat
                   position within +/-6 months of each dune-line survey date,
                   end minus start, domain mean. The like-for-like quantity
                   for two surveys.
    Both read from the stored products, not computed here. Shown as the net
    change over the survey interval (m/yr) so the four windows share one axis;
    the metres are in domain_comparison.csv.

    Until 09-18 the shoreline was ALSO drawn as the CoastSat LRR, an OLS
    through ~250 dates, and the correlations reported against both. That is
    not a two-survey quantity; it stays the model's scoring target in
    3-rates/coastsat/lrr/ and in model_vs_observed/vs_shoreline/. The helpers below
    (window_mean, endpoint_by_transect, KNOWN_SURVEY_DATES) are kept: the
    stored CoastSat endpoint product is built with them.

SURVEY DATES
    The dune-line files carry no date. KNOWN_SURVEY_DATES holds them by line
    VINTAGE: 1984-09-19 and 1997-10-12 from the Henderson USGS metadata on
    D:, 2004-05-25 and 2009-05-30 from the Google Earth captures (Hannah,
    2026-09-15), and None for the 2023 NOAA set until its flight date is
    known. A None is centred mid-year of the line's year, PROVENANCE.md says
    so, and a sensitivity block reports how far the endpoint rate moves for
    a +/- 6 month shift of that centre. A period year reaches its vintage
    through hat_topo_version.DUNE_LINE_FOR_YEAR (2010 reads the 2009 line,
    2024 the 2023 one).

METHOD
    Every raw dune file is built by duneline_to_raw_offsets.py since
    2026-09-15 (1984 and 2004 were ArcGIS exports until that afternoon, one
    metre landward of the exact crossing; see raw_offsets/PROVENANCE.md), so
    a change between any two years carries no method term.

OUTPUT   data/hatteras_init/5-scr/4-comparisons/shoreline_vs_duneline/endpoint_net_change/<start>_<end>/
         (was 4-comparisons/coastsat_vs_duneline/ until 2026-09-19; the
         alongshore figures are in METRES with the beach-width gap since then)
             scatter_dune_vs_coastsat.png     shoreline vs dune, 1:1
             alongshore_dune_vs_coastsat.png  the two net changes by domain
             supporting/                      the PDFs, CAPTIONS.md,
                 domain_comparison.csv            one row per GIS domain (m and m/yr)
                 PROVENANCE.md
         data/hatteras_init/5-scr/4-comparisons/shoreline_vs_duneline/endpoint_net_change/
             alongshore_four_windows.png      every window on one y axis,
                                              stacked full width (--grid;
                                              --layout grid for the 2 x 2)

USAGE
    python coastsat_vs_duneline.py --start-year 1984 --end-year 2004
        # the dates come from the stored products
    python coastsat_vs_duneline.py --grid                 # every window, stacked
    python coastsat_vs_duneline.py --grid --layout grid   # the 2 x 2 by period
==============================================================================
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = next(_p for _p in Path(__file__).resolve().parents
                    if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from site_layer.hat_observed_rates import (COASTSAT_ENDPOINT_VS_DUNELINE,  # noqa: E402
                                          COASTSAT_TIMESERIES,  # noqa: E402
                                SCR_ROOT, coastsat_endpoint_csv, dune_endpoint_csv,
                                transect_lookup)
from site_layer.hat_topo_version import dune_line_for_year, dune_raw_file_for_year  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

from site_layer.hat_figure_style import (  # noqa: E402
    compare_header,C, C_1984, C_1997, DOMAIN_AXIS_LABEL,  # noqa: E402
                              INK_MUTED, _title, apply_style, caption,
                              figsize, open_frame, save, structures,
                              support_dir, town_bands)

OUT_ROOT = COASTSAT_ENDPOINT_VS_DUNELINE
GIS_FIRST, GIS_LAST = 1, 90
DAYS_PER_YEAR = 365.25
SIX_MONTHS_DAYS = 182.625

# Colours (Hannah, 2026-09-15, third pass): the two CoastSat estimators are
# the SAME feature measured two ways, so they share the blue family, one line
# each -- the LRR dark (the RdBu blue pole the house uses for the sea side),
# the endpoint a lighter blue and dashed so the pair still separates in
# greyscale. The dune line is the house RdBu red against them, the pair the
# other figures already use, so the red/blue contrast reads at a glance
# (Hannah, 2026-09-15). Here the pair means FEATURE, dune vs shoreline, not
# vintage and not sign; the caption says so. Grey, purple, sand brown, orange
# and REF green were tried and rejected, the green as too dark to see.
C_LRR = C_1997                  # "#2166ac"
C_ENDPOINT = "#74a9cf"          # PuBu mid blue (unused since 09-18)
# Since 2026-09-18 the shoreline is ONE line, its net change at the dune
# dates, drawn in the dark blue the LRR had.
C_SHORE = C_LRR
C_DUNE = C_1984                 # "#b2182b"

# Keyed by LINE VINTAGE (the year in the geojson name), not by period year: a
# period finds its vintage through hat_topo_version.DUNE_LINE_FOR_YEAR.
# 1984, 1997: the Henderson USGS metadata on D: (Calendar_Date). 2004, 2009:
# Google Earth capture dates (Hannah, 2026-09-15); the raw_GE frames carry
# none. 2023: NOAA NGS imagery under D:\Hatteras_GIS\Aerial\2023 whose
# metadata gives only the 2015-2023 series extent, so None until Hannah
# supplies the flight date; a None is centred mid-year and flagged.
KNOWN_SURVEY_DATES = {1984: "1984-09-19", 1997: "1997-10-12",
                      2004: "2004-05-25", 2009: "2009-05-30", 2023: None}


# -----------------------------------------------------------------------------
# dune side
# -----------------------------------------------------------------------------
def dune_position_by_domain(year: int) -> pd.Series:
    """Mean ORIG_LEN per GIS domain, first row per transect, as the hindcast
    loader (`load_absolute_dune_distance`) reads it. Grows LANDWARD."""
    path = dune_raw_file_for_year(year)     # the vintage that stands for `year`
    raw = pd.read_csv(path, encoding="utf-8-sig")
    per_transect = raw.drop_duplicates(subset=["domain_id", "LineID"])
    means = per_transect.groupby("domain_id")["ORIG_LEN"].mean()
    return means.reindex(range(GIS_FIRST, GIS_LAST + 1))


# -----------------------------------------------------------------------------
# coastsat side
# -----------------------------------------------------------------------------
def load_chainage(transect_id: str):
    site = transect_id.rsplit("_", 1)[0]
    path = COASTSAT_TIMESERIES / f"{site}_timeseries" / f"{transect_id}.csv"
    if not path.is_file():
        return None
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    df = df.rename(columns={"dates UTC": "date", "chainage (m)": "chainage"})
    df["date"] = pd.to_datetime(df["date"], utc=True)
    return df.dropna(subset=["chainage"])


def window_mean(df: pd.DataFrame, centre: datetime, half_days: float):
    lo = centre - timedelta(days=half_days)
    hi = centre + timedelta(days=half_days)
    sel = df[(df["date"] >= lo) & (df["date"] <= hi)]
    if sel.empty:
        return np.nan, 0, np.nan, pd.NaT, pd.NaT
    ch = sel["chainage"]
    return (float(ch.mean()), int(ch.size), float(ch.std(ddof=0)),
            sel["date"].min(), sel["date"].max())


def endpoint_by_transect(lookup: pd.DataFrame, d0: datetime, d1: datetime,
                         half_days: float, cache: dict) -> pd.DataFrame:
    years = (d1 - d0).days / DAYS_PER_YEAR
    rows = []
    for tid, dom in zip(lookup["transect_id"], lookup["domain_number"]):
        if tid not in cache:
            cache[tid] = load_chainage(tid)
        df = cache[tid]
        if df is None:
            continue
        m0, n0, s0, f0, l0 = window_mean(df, d0, half_days)
        m1, n1, s1, f1, l1 = window_mean(df, d1, half_days)
        rows.append(dict(transect_id=tid, domain_number=int(dom),
                         mean_start_m=m0, n_start=n0, sd_start_m=s0,
                         first_obs_start=f0, last_obs_start=l0,
                         mean_end_m=m1, n_end=n1, sd_end_m=s1,
                         first_obs_end=f1, last_obs_end=l1,
                         change_m=m1 - m0, endpoint_rate_m_yr=(m1 - m0) / years))
    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# figures
# -----------------------------------------------------------------------------
def _fit_stats(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if x.size < 3:
        return dict(n=int(x.size), r=np.nan, slope=np.nan, intercept=np.nan,
                    rmse=np.nan, bias=np.nan)
    slope, intercept = np.polyfit(x, y, 1)
    r = np.corrcoef(x, y)[0, 1]
    return dict(n=int(x.size), r=float(r), slope=float(slope),
                intercept=float(intercept),
                rmse=float(np.sqrt(np.mean((y - x) ** 2))),
                bias=float(np.mean(y - x)))


def scatter_figure(dom: pd.DataFrame, out: Path, start: int, end: int,
                   st: dict, v0: int, v1: int) -> None:
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=figsize("single", aspect=1.0),
                           constrained_layout=True)
    x = dom["dune_rate_m_yr"].to_numpy(dtype=float)
    y = dom["cs_endpoint_m_yr"].to_numpy(dtype=float)
    lim = np.nanmax(np.abs(np.concatenate([x, y]))) * 1.05
    ax.plot([-lim, lim], [-lim, lim], color=INK_MUTED, lw=0.6, ls="--", zorder=1)
    ax.axhline(0, color=C["GRID"], lw=0.6, zorder=0)
    ax.axvline(0, color=C["GRID"], lw=0.6, zorder=0)
    ax.scatter(x, y, s=12, color=C["BASE"], edgecolor="none", zorder=3)
    resid = y - x
    for j in np.argsort(-np.abs(np.nan_to_num(resid)))[:6]:
        if np.isfinite(resid[j]):
            ax.annotate(str(int(dom["gis"].iloc[j])), (x[j], y[j]),
                        xytext=(3, 3), textcoords="offset points",
                        fontsize=6.5, color=INK_MUTED)
    if np.isfinite(st["slope"]):
        xx = np.array([-lim, lim])
        ax.plot(xx, st["slope"] * xx + st["intercept"], color=C["ACCENT"],
                lw=0.9, zorder=2)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.set_xlabel("Dune line, net change rate (m/yr)")
    ax.set_ylabel("CoastSat shoreline, net change rate (m/yr)")
    caption(fig, (
        f"Per-domain net change of the digitized dune line against the CoastSat "
        f"shoreline over the same interval, {v0}–{v1} (standing in for "
        f"{start}–{end}), n = {st['n']} domains, both as net change over the "
        "survey interval, seaward positive. The shoreline is the mean CoastSat "
        "position within six months of each dune-line image date, differenced. "
        f"Dashed: 1:1; purple: least-squares fit (r = {st['r']:.2f}, slope = "
        f"{st['slope']:.2f}, RMSE = {st['rmse']:.2f} m/yr, bias shoreline − dune "
        f"= {st['bias']:+.2f} m/yr). The six domains farthest from 1:1 are "
        "labelled. See PROVENANCE.md for the survey dates."))
    # The window is IN the stem (Hannah, 2026-09-21): five windows wrote this
    # same basename, so a figure lifted out of its folder could not be told
    # from the other four.
    save(fig, out / f"coastsat_endpoint_vs_duneline_{start}_{end}_scatter", close=True)


# The alongshore figures are in METRES since 2026-09-19 (Hannah: one form for
# every shoreline-vs-dune figure, net change with the beach-width gap).
Y_TICK = 20.0
Y_LABEL = "Net change in position (m)"
# The gap between the two lines: widened solid grey, narrowed hatched. Shared
# by net_change_vs_duneline.py and total_change_vs_duneline.py.
C_GAP = "0.86"
C_NARROW = "0.55"
NARROW_HATCH = "////"
STRUCTURE_LABEL_PT_GRID = 4.0   # the 2 x 2, whose panels are half the width


def _half(*arrays) -> float:
    """Symmetric y limit: the largest |value| plus 5 m, up to the next 10 m."""
    v = np.concatenate([np.asarray(a, dtype=float).ravel() for a in arrays])
    return float(np.ceil((np.nanmax(np.abs(v)) + 5.0) / 10.0) * 10.0)


def shade_beach_width(ax, x, shore, dune, zorder=3):
    """The space between the shoreline and dune-line changes: solid grey where
    the beach WIDENED (shoreline change > dune-line change), hatched where it
    narrowed."""
    x, shore, dune = (np.asarray(a, dtype=float) for a in (x, shore, dune))
    width = shore - dune
    ax.fill_between(x, dune, shore, where=width >= 0, interpolate=True,
                    color=C_GAP, lw=0, zorder=zorder)
    ax.fill_between(x, dune, shore, where=width < 0, interpolate=True,
                    facecolor="white", edgecolor=C_NARROW, hatch=NARROW_HATCH,
                    lw=0, zorder=zorder)


def beach_width_handles():
    from matplotlib.patches import Patch
    return [Patch(facecolor=C_GAP, lw=0, label="Beach widened"),
            Patch(facecolor="white", edgecolor=C_NARROW, hatch=NARROW_HATCH, lw=0,
                  label="Beach narrowed")]


def _legend_handles():
    return [
        Line2D([], [], color=C_SHORE, lw=1.1, label="Shoreline change (CoastSat endpoint)"),
        Line2D([], [], color=C_DUNE, lw=1.1, label="Total dune line change (measured)"),
    ] + beach_width_handles()


def draw_alongshore(ax, dom: pd.DataFrame, half: float, label: bool = True,
                    label_pt: float = 6.5, window=None) -> None:
    """One alongshore panel in the house form: village bands, groin and piers,
    open frame, y grid, symmetric limits. TWO lines since 2026-09-18, both net
    change over the same survey interval: the CoastSat shoreline blue, the
    dune line red; in METRES with the beach-width gap shaded since 2026-09-19.
    Call it after the legend is placed (structures() tests its labels against
    the layout)."""
    from matplotlib.ticker import MultipleLocator

    g = dom["gis"].to_numpy(dtype=float)
    ax.set_xlim(GIS_FIRST - 0.5, GIS_LAST + 0.5)
    ax.set_ylim(-half, half)
    town_bands(ax, label=label)
    ax.axhline(0, color=INK_MUTED, lw=0.6, zorder=2)
    shore = dom["cs_endpoint_change_m"].to_numpy(dtype=float)
    dune = dom["dune_change_m"].to_numpy(dtype=float)
    shade_beach_width(ax, g, shore, dune)
    ax.plot(g, dune, color=C_DUNE, lw=1.1, zorder=5)
    ax.plot(g, shore, color=C_SHORE, lw=1.1, zorder=5)
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_minor_locator(MultipleLocator(5))
    ax.yaxis.set_major_locator(MultipleLocator(Y_TICK))
    ax.yaxis.grid(True, zorder=0)
    ax.set_axisbelow(True)
    open_frame(ax)
    structures(ax, label=label, label_pt=label_pt)
    # the offshore shoals, and the fills placed inside the window, as the
    # 3-rates figures mark them (2026-09-19)
    sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "input_prep"
                           / "5-scr" / "lib"))
    import scr_paths  # noqa: F401  (5-scr sibling modules onto sys.path)
    import coastsat_lrr_windows as cw
    cw.draw_shoals(ax, label=label, label_pt=label_pt)
    if window is not None:
        fills = cw.fills_in(*window)
        if fills:
            cw.draw_fills(ax, fills, half, label_pt=label_pt)


def _caption_body() -> str:
    return ("Both lines are NET CHANGE in metres over the same survey interval, "
            "seaward positive. Red: the digitized dune line, end line minus start "
            "line. Blue: the CoastSat shoreline, the mean satellite position within "
            "six months of each dune-line image date, end minus start. Both are read "
            "from the stored products in 5-scr/3-rates (duneline/endpoint, "
            "coastsat/endpoint). The space between them is beach-width change "
            "(shoreline minus dune line): solid grey where the beach widened, "
            "hatched where it narrowed. "
            "Red and blue here mark the two features, not the sign or the vintage. "
            "Hatched amber boxes mark the offshore shoals and black bars above the "
            "panel the beach fills placed in the window. "
            "Bands mark Buxton, Avon and the Tri-Village; the solid hairline is the "
            "Buxton groin, the dotted ones the Avon and Rodanthe piers. Domain 1 is "
            "Cape Point, 90 is Pea Island.")


def alongshore_figure(dom: pd.DataFrame, out: Path, start: int, end: int,
                      v0: int, v1: int) -> None:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=figsize("double", aspect=0.38),
                           constrained_layout=True)
    half = _half(dom["cs_endpoint_change_m"], dom["dune_change_m"])
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    ax.set_ylabel(Y_LABEL)
    fig.legend(handles=_legend_handles(), loc="outside lower center", ncol=4,
               frameon=False)
    draw_alongshore(ax, dom, half, window=(start, end))
    caption(fig, f"Net change of the dune line and the CoastSat shoreline, "
                 f"{v0}–{v1} (standing in for {start}–{end}), by GIS domain. "
                 + _caption_body())
    # Both sides are OBSERVED here -- two snapshots differenced, no rate
    # anywhere -- which is what separates this folder from total_change/.
    compare_header(fig, [
        f"{start}–{end}   ·   shoreline: CoastSat, mean position ±6 months about "
        "each dune-line date",
        f"dune line: {v0} → {v1} lines,  measured"])
    save(fig, out / f"coastsat_endpoint_vs_duneline_{start}_{end}_alongshore", close=True)


def _chains(windows):
    """Windows linked end-to-start: [(1984,2004),(2004,2024)],
    [(1996,2010),(2010,2024)] -- the same rule as coastsat_lrr_windows.py."""
    rest = sorted(windows)
    chains = []
    while rest:
        chain = [rest.pop(0)]
        while True:
            nxt = next((w for w in rest if w[0] == chain[-1][1]), None)
            if nxt is None:
                break
            rest.remove(nxt)
            chain.append(nxt)
        chains.append(chain)
    return chains


def four_windows_figure(layout: str = "column") -> Path:
    """Every model window on ONE y axis, one full-width panel per window in
    chain order (`layout="column"`), or the 2 x 2 by period (`"grid"`). Reads
    each window's supporting/domain_comparison.csv; run the windows first.
    The context window 1996_2024 is left out: it is not a model window."""
    import matplotlib.pyplot as plt

    windows, frames = [], {}
    for d in sorted(OUT_ROOT.iterdir()):
        s_, _, e_ = d.name.partition("_")
        f = d / "supporting" / "domain_comparison.csv"
        if s_.isdigit() and e_.isdigit() and f.is_file() and (int(s_), int(e_)) != (1996, 2024):
            windows.append((int(s_), int(e_)))
            frames[(int(s_), int(e_))] = pd.read_csv(f)
    if not windows:
        sys.exit(f"no windows under {OUT_ROOT}; run the comparison first")
    half = _half(*[frames[w][c] for w in windows
                   for c in ("cs_endpoint_change_m", "dune_change_m")])

    chains = _chains(windows)
    grid = layout == "grid" and len(chains) == 2 and all(len(c) == 2 for c in chains)
    if grid:
        nrow, ncol = 2, 2
        cells = [(r, c, chain[r]) for r in range(2) for c, chain in enumerate(chains)]
        fig, axes = plt.subplots(nrow, ncol, sharex=True, sharey=True,
                                 figsize=figsize("double", height=4.9),
                                 constrained_layout=True)
        label_pt = STRUCTURE_LABEL_PT_GRID
    else:
        ordered = [w for chain in chains for w in chain]
        nrow, ncol = len(ordered), 1
        cells = [(i, 0, w) for i, w in enumerate(ordered)]
        fig, axes = plt.subplots(nrow, ncol, sharex=True, sharey=True,
                                 figsize=figsize("double", height=min(2.0 * nrow + 0.9, 9.4)),
                                 constrained_layout=True, squeeze=False)
        label_pt = 6.5
    axes = np.asarray(axes).reshape(nrow, ncol)
    for ax in axes[-1, :]:
        ax.set_xlabel(DOMAIN_AXIS_LABEL)
    fig.supylabel(Y_LABEL, fontsize=9)
    fig.legend(handles=_legend_handles(), loc="outside lower center", ncol=4,
               frameon=False)
    for i, (r, c, (start, end)) in enumerate(cells):
        ax = axes[r, c]
        draw_alongshore(ax, frames[(start, end)], half, label=(i == 0),
                        label_pt=label_pt, window=(start, end))
        _title(ax, i, f"{start}–{end}")
        if c > 0:
            ax.tick_params(labelleft=False)
    wins = ", ".join(f"{a}–{b}" for _, _, (a, b) in cells)
    layout = ("a 2 x 2 with the 1984-start period in the left column and the "
              "1996-start period in the right, the earlier window above"
              if grid else "one full-width panel per window, the 1984-start "
              "pair above the 1996-start pair")
    caption(fig, (f"Net change of the dune line and the CoastSat shoreline by GIS "
                  f"domain for the {len(cells)} hindcast windows ({wins}), "
                  f"{layout}, all on one y axis (±{half:g} m, the largest "
                  f"value over every window plus 5 m). " + _caption_body()
                  + " Each window's survey dates and statistics are in its "
                  "own supporting/PROVENANCE.md."))
    return save(fig, OUT_ROOT / "coastsat_endpoint_vs_duneline_four_windows_alongshore",
                close=True)[0]


# -----------------------------------------------------------------------------
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="net dune-line change against net CoastSat shoreline change")
    ap.add_argument("--start-year", type=int, default=1984)
    ap.add_argument("--end-year", type=int, default=2004)
    ap.add_argument("--half-window-days", type=float, default=SIX_MONTHS_DAYS,
                    help="half-width of the CoastSat window, for the assumed-"
                         "date sensitivity only (the product fixes the main one)")
    ap.add_argument("--grid", action="store_true",
                    help="draw every model window on disk as one figure "
                         "(alongshore_four_windows) and exit; no window is run")
    ap.add_argument("--layout", choices=("column", "grid"), default="column",
                    help="with --grid: one full-width panel per window "
                         "(default), or the 2 x 2 by period")
    a = ap.parse_args(argv)
    if a.grid:
        print(f"-> {four_windows_figure(a.layout)}")
        return 0

    s, e = a.start_year, a.end_year
    # BOTH SIDES FROM THE STORED NET-CHANGE PRODUCTS (2026-09-18)
    du_dom = pd.read_csv(dune_endpoint_csv(s, e, "domain")).set_index("domain_number")
    du_tr = pd.read_csv(dune_endpoint_csv(s, e, "transect"))
    cs_dom = pd.read_csv(coastsat_endpoint_csv(s, e, "domain")).set_index("domain_number")
    cs_tr = pd.read_csv(coastsat_endpoint_csv(s, e, "transect"))
    meta = du_tr.iloc[0]
    v0, v1 = int(meta["start_vintage"]), int(meta["end_vintage"])
    d0 = datetime.fromisoformat(meta["start_date"]).replace(tzinfo=timezone.utc)
    d1 = datetime.fromisoformat(meta["end_date"]).replace(tzinfo=timezone.utc)
    assumed = {k for k, c in (("start", "start_date_assumed"), ("end", "end_date_assumed"))
               if bool(meta[c])}
    years = float(meta["interval_yr"])
    print(f"survey dates  {d0.date()} -> {d1.date()}  ({years:.2f} yr)"
          + ("   [END DATE ASSUMED]" if "end" in assumed else ""))

    apply_style()
    out = OUT_ROOT / f"{s}_{e}"
    out.mkdir(parents=True, exist_ok=True)
    sup = support_dir(out)
    old = sup / "transect_coastsat_endpoint.csv"   # now the stored product
    if old.is_file():
        old.unlink()

    pos = du_tr.groupby("domain_number")[[f"position_{v0}_m", f"position_{v1}_m"]].mean()
    dom = pd.DataFrame({"gis": range(GIS_FIRST, GIS_LAST + 1)}).set_index("gis")
    dom["dune_start_m"] = pos[f"position_{v0}_m"]
    dom["dune_end_m"] = pos[f"position_{v1}_m"]
    dom["dune_change_m"] = du_dom["mean_change_m"]
    dom["dune_rate_m_yr"] = du_dom["mean_rate_m_yr"]
    dom["cs_endpoint_change_m"] = cs_dom["mean_change_m"]
    dom["cs_endpoint_m_yr"] = cs_dom["mean_rate_m_yr"]
    dom["n_obs_start"] = cs_dom["median_n_start"]
    dom["n_obs_end"] = cs_dom["median_n_end"]
    dom["n_transects"] = cs_dom["n_transects"]
    dom["beach_width_change_m"] = dom["cs_endpoint_change_m"] - dom["dune_change_m"]
    dom["dune_minus_endpoint_m_yr"] = dom["dune_rate_m_yr"] - dom["cs_endpoint_m_yr"]
    dom = dom.reset_index()
    dom.to_csv(sup / "domain_comparison.csv", index=False, float_format="%.3f")

    st = _fit_stats(dom["dune_rate_m_yr"].to_numpy(), dom["cs_endpoint_m_yr"].to_numpy())

    span = {}
    for key, d in (("start", d0), ("end", d1)):
        f = pd.to_datetime(cs_tr[f"first_obs_{key}"]).median()
        l = pd.to_datetime(cs_tr[f"last_obs_{key}"]).median()
        lo = (d - timedelta(days=SIX_MONTHS_DAYS)).date()
        hi = (d + timedelta(days=SIX_MONTHS_DAYS)).date()
        span[key] = (lo, hi, f.date(), l.date(),
                     (f.date() - lo).days > 45 or (hi - l.date()).days > 45)

    # sensitivity of the CoastSat net change to an assumed survey date
    sens = []
    if assumed:
        lookup = pd.read_csv(transect_lookup())
        lookup = lookup[lookup["domain_number"].between(GIS_FIRST, GIS_LAST)]
        cache: dict = {}
        for label, shift in (("-6 mo", -SIX_MONTHS_DAYS), ("0", 0.0),
                             ("+6 mo", SIX_MONTHS_DAYS)):
            dd0 = d0 + timedelta(days=shift if "start" in assumed else 0)
            dd1 = d1 + timedelta(days=shift if "end" in assumed else 0)
            ep = endpoint_by_transect(lookup, dd0, dd1, a.half_window_days, cache)
            e_dom = ep.groupby("domain_number")["endpoint_rate_m_yr"].mean()
            ss = _fit_stats(
                dom.set_index("gis")["dune_rate_m_yr"].reindex(e_dom.index).to_numpy(),
                e_dom.to_numpy())
            sens.append((label, dd0.date(), dd1.date(), float(e_dom.mean()), ss))

    scatter_figure(dom, out, s, e, st, v0, v1)
    alongshore_figure(dom, out, s, e, v0, v1)
    isl = dom[["dune_rate_m_yr", "cs_endpoint_m_yr", "dune_change_m",
               "cs_endpoint_change_m"]].mean()
    date_rows = []
    for key, v, d in (("start", v0, d0), ("end", v1, d1)):
        yr = s if key == "start" else e
        stand_in = f" — the {v} line standing in for {yr}" if v != yr else ""
        src = ("**ASSUMED 1 July**; no flight date known for this line"
               if key in assumed else f"`coastsat_vs_duneline.KNOWN_SURVEY_DATES`")
        date_rows.append(f"| {yr} | {d.date()} | {src}{stand_in} |")
    lines = [
        f"# Dune line vs CoastSat shoreline, {s}-{e} (net change)",
        "",
        f"Written {datetime.now():%Y-%m-%d %H:%M} by "
        "`scripts/input_prep/5-scr/4-comparisons/shoreline_vs_duneline/coastsat_vs_duneline.py`.",
        "",
        "**Both sides are NET CHANGE between the same two dates** (2026-09-18, "
        "Hannah: the comparison is net position change on both sides). The "
        "CoastSat LRR, which this folder also drew until then, is not a like-for-"
        "like quantity for two surveys; it stays the model's scoring target in "
        "`3-rates/coastsat/lrr/`.",
        "",
        "## Inputs (the stored products)",
        "",
        f"* dune line: `3-rates/duneline/endpoint/{s}_{e}/` (the {v0} and {v1} lines).",
        f"* CoastSat shoreline: `3-rates/coastsat/endpoint/{s}_{e}/` (mean position "
        "within ±6 months of each dune-line date, per transect, domain mean).",
        "",
        "## Survey dates",
        "",
        "| line | date | source |",
        "|---|---|---|",
        *date_rows,
        "",
        f"Survey interval {years:.2f} yr. Seaward positive in every column.",
        "",
        "## Island-wide means",
        "",
        "| | net change (m) | as a rate (m/yr) |",
        "|---|---|---|",
        f"| dune line | {isl['dune_change_m']:+.1f} | {isl['dune_rate_m_yr']:+.2f} |",
        f"| CoastSat shoreline | {isl['cs_endpoint_change_m']:+.1f} | {isl['cs_endpoint_m_yr']:+.2f} |",
        "",
        "## Agreement, per domain (shoreline y against dune x, m/yr)",
        "",
        "| n | r | slope | intercept | RMSE | bias (y − x) |",
        "|---|---|---|---|---|---|",
        f"| {st['n']} | {st['r']:.2f} | {st['slope']:.2f} | {st['intercept']:.2f} | "
        f"{st['rmse']:.2f} | {st['bias']:+.2f} |",
        "",
        "## Window occupancy (CoastSat)",
        "",
        f"Median positions per transect inside the start window "
        f"{int(cs_tr['n_start'].median())}, the end window {int(cs_tr['n_end'].median())}; "
        f"empty windows: start {int((cs_tr['n_start'] == 0).sum())}, "
        f"end {int((cs_tr['n_end'] == 0).sum())}, of {len(cs_tr)} transects.",
        "",
        "| window | asked for | median first obs | median last obs | truncated |",
        "|---|---|---|---|---|",
    ]
    for key in ("start", "end"):
        lo, hi, f, l, trunc = span[key]
        lines.append(f"| {key} | {lo} to {hi} | {f} | {l} | {'**yes**' if trunc else 'no'} |")
    if any(v[4] for v in span.values()):
        lines += ["", "A truncated window does not average a full seasonal cycle. "
                  "The CoastSat record begins 1984-09-21 on most transects, so a "
                  "window centred on the 1984-09-19 photo holds only the autumn and "
                  "winter after it."]
    if sens:
        lines += ["", "## Sensitivity of the CoastSat net change to the assumed date", "",
                  "| centre shift | start | end | island mean (m/yr) | r vs dune | slope | RMSE |",
                  "|---|---|---|---|---|---|---|"]
        for label, s0, s1, m, ss in sens:
            lines.append(f"| {label} | {s0} | {s1} | {m:.2f} | {ss['r']:.2f} | "
                         f"{ss['slope']:.2f} | {ss['rmse']:.2f} |")
    lines += ["", "## Read this before quoting it", "",
              "* The dune line and the waterline are different features; a gap "
              "between them is beach-width change as much as disagreement.",
              "* One storm inside a ±6-month window moves that end."]
    (sup / "PROVENANCE.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"island means  dune {isl['dune_rate_m_yr']:+.2f}  coastsat "
          f"{isl['cs_endpoint_m_yr']:+.2f} m/yr   (net {isl['dune_change_m']:+.1f} / "
          f"{isl['cs_endpoint_change_m']:+.1f} m)")
    print(f"  shoreline vs dune  n={st['n']} r={st['r']:.2f} slope={st['slope']:.2f} "
          f"rmse={st['rmse']:.2f} bias={st['bias']:+.2f}")
    for label, s0, s1, m, ss in sens:
        print(f"  sens {label:6s} {s0} {s1}  mean={m:+.2f}  r={ss['r']:.2f}")
    print(f"-> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
