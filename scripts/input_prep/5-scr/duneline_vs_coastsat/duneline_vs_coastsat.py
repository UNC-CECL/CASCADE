r"""
duneline_vs_coastsat.py
==============================================================================
Does the digitized dune line move with the CoastSat shoreline?

WHAT IS COMPARED
    Dune side      two digitized dune lines (start year, end year), read from
                   2-brie-offset/raw_offsets/<year>_duneline_offset_raw.csv as
                   the hindcast loader reads them: first row per transect,
                   mean over the ~5 transects of each GIS domain. ORIG_LEN is
                   the station from the OFFSHORE datum line, so it grows
                   landward; the change is negated so that, like CoastSat,
                   seaward is positive and a negative rate is retreat.
    Shoreline side the CoastSat waterline on the same 90 domains, through
                   transect_domain_lookup.csv, as TWO estimators:
                     lrr       the per-transect OLS slope already on disk in
                               coastsat_lrr/<start>_<end>/transect_lrr_full.csv,
                               the quantity the model is graded against;
                     endpoint  the mean chainage inside +/- HALF_WINDOW days of
                               each SURVEY date, differenced and divided by
                               the survey interval. This is the like-for-like
                               quantity for two surveys. A symmetric one-year
                               window averages one full seasonal cycle.
                   Both are averaged to the domain (mean over transects, as
                   domain_lrr_summary.csv does).

WHY BOTH  (Hannah, 2026-09-15)
    A dune line is two moments; an OLS slope through ~250 satellite dates is
    not the same quantity, so the endpoint rate is the fair comparison and the
    LRR is the one the rest of 5-scr uses. Report both, and the gap.

SURVEY DATES
    The dune-line files carry no date. 1984 is the 1984-09-19 USGS photo
    (Henderson, D:\Hatteras_GIS\Aerial\1984_henderson\1984_metadata). The 2004
    line was traced from Google Earth captures whose file carries no date; the
    capture is 2004-05-25 (Hannah, 2026-09-15). Both are in KNOWN_SURVEY_DATES.
    A year with no known date is centred mid-year, PROVENANCE.md says so, and
    a sensitivity block reports how far the endpoint rate moves for a
    +/- 6 month shift of that centre.

WHAT CANCELS AND WHAT DOES NOT
    Both the 1984 and 2004 raw files are ArcGIS exports, so the 1 m GIS-vs-
    shapely convention (raw_offsets/PROVENANCE.md) cancels here. It will NOT
    cancel against a shapely-built 2024 file; that leg needs both lines
    through duneline_to_raw_offsets.py or a stated 1 m correction.

OUTPUT   data/hatteras_init/5-scr/duneline_vs_coastsat/<start>_<end>/
             scatter_dune_vs_coastsat.png     a. vs LRR  b. vs endpoint
             alongshore_dune_vs_coastsat.png  the three rates by domain
             supporting/                      the PDFs, CAPTIONS.md,
                 domain_comparison.csv            one row per GIS domain
                 transect_coastsat_endpoint.csv   the window means per transect
                 PROVENANCE.md

USAGE
    python duneline_vs_coastsat.py --start-year 1984 --end-year 2004 \
        # dates come from KNOWN_SURVEY_DATES; --start-date/--end-date override
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

from hat_observed_rates import (COASTSAT_TIMESERIES, SCR_ROOT,  # noqa: E402
                                lrr_csv, transect_lookup)
from matplotlib.lines import Line2D  # noqa: E402

from hat_figure_style import (C, C_1984, C_1997, DOMAIN_AXIS_LABEL,  # noqa: E402
                              INK_MUTED, _title, apply_style, caption,
                              figsize, open_frame, save, structures,
                              support_dir, town_bands)

RAW_DIR = PROJECT_ROOT / "data" / "hatteras_init" / "2-brie-offset" / "raw_offsets"
OUT_ROOT = SCR_ROOT / "duneline_vs_coastsat"
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
C_ENDPOINT = "#74a9cf"          # PuBu mid blue
C_DUNE = C_1984                 # "#b2182b"

# 1984: the Henderson metadata. 2004: the Google Earth capture date (Hannah,
# 2026-09-15); the raw_GE frames themselves carry none.
KNOWN_SURVEY_DATES = {1984: "1984-09-19", 2004: "2004-05-25"}


# -----------------------------------------------------------------------------
# dune side
# -----------------------------------------------------------------------------
def dune_position_by_domain(year: int) -> pd.Series:
    """Mean ORIG_LEN per GIS domain, first row per transect, as the hindcast
    loader (`load_absolute_dune_distance`) reads it. Grows LANDWARD."""
    path = RAW_DIR / f"{year}_duneline_offset_raw.csv"
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
                   stats: dict) -> None:
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=figsize("double", aspect=0.5),
                             constrained_layout=True, sharex=True, sharey=True)
    pairs = [("cs_lrr_m_yr", "CoastSat LRR (m/yr)"),
             ("cs_endpoint_m_yr", "CoastSat endpoint rate (m/yr)")]
    x = dom["dune_rate_m_yr"].to_numpy()
    allv = np.concatenate([x] + [dom[c].to_numpy() for c, _ in pairs])
    lim = np.nanmax(np.abs(allv)) * 1.05
    for i, (ax, (col, lab)) in enumerate(zip(axes, pairs)):
        y = dom[col].to_numpy()
        ax.plot([-lim, lim], [-lim, lim], color=INK_MUTED, lw=0.6, ls="--",
                zorder=1)
        ax.axhline(0, color=C["GRID"], lw=0.6, zorder=0)
        ax.axvline(0, color=C["GRID"], lw=0.6, zorder=0)
        ax.scatter(x, y, s=12, color=C["BASE"], edgecolor="none", zorder=3)
        # label the far-from-1:1 domains so a reader can find them
        resid = y - x
        far = np.argsort(-np.abs(np.nan_to_num(resid)))[:6]
        for j in far:
            if np.isfinite(resid[j]):
                ax.annotate(str(int(dom["gis"].iloc[j])), (x[j], y[j]),
                            xytext=(3, 3), textcoords="offset points",
                            fontsize=6.5, color=INK_MUTED)
        st = stats[col]
        if np.isfinite(st["slope"]):
            xx = np.array([-lim, lim])
            ax.plot(xx, st["slope"] * xx + st["intercept"], color=C["ACCENT"],
                    lw=0.9, zorder=2)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.set_xlabel("Dune-line rate (m/yr)")
        ax.set_ylabel(lab)
        _title(ax, i, lab.replace(" (m/yr)", ""))
    caption(fig, (
        f"Per-domain rate of the digitized dune line, {start}-{end}, against "
        f"the CoastSat shoreline on the same GIS domains (n = {stats['cs_lrr_m_yr']['n']}). "
        f"Seaward positive. Dashed: 1:1; purple: least-squares fit. "
        f"(a) the per-transect OLS slope over the window, averaged per domain "
        f"(r = {stats['cs_lrr_m_yr']['r']:.2f}, slope = {stats['cs_lrr_m_yr']['slope']:.2f}, "
        f"RMSE = {stats['cs_lrr_m_yr']['rmse']:.2f} m/yr). "
        f"(b) the endpoint rate from the mean CoastSat position in a one-year "
        f"window centred on each survey date "
        f"(r = {stats['cs_endpoint_m_yr']['r']:.2f}, slope = {stats['cs_endpoint_m_yr']['slope']:.2f}, "
        f"RMSE = {stats['cs_endpoint_m_yr']['rmse']:.2f} m/yr). "
        f"The six domains farthest from 1:1 are labelled. See PROVENANCE.md "
        f"for the survey dates and what the {end} date assumes."))
    save(fig, out / "scatter_dune_vs_coastsat", close=True)


def alongshore_figure(dom: pd.DataFrame, out: Path, start: int, end: int) -> None:
    """The house alongshore panel (coastsat_lrr_windows.py): full-height
    village bands, the groin and piers named along their lines, open frame,
    y grid, symmetric y limits on the 2 m tick. Three lines: the two CoastSat
    estimators in one colour family, the dune line warm against them."""
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator

    fig, ax = plt.subplots(figsize=figsize("double", aspect=0.38),
                           constrained_layout=True)
    g = dom["gis"].to_numpy(dtype=float)
    lrr = dom["cs_lrr_m_yr"].to_numpy(dtype=float)
    ept = dom["cs_endpoint_m_yr"].to_numpy(dtype=float)
    dune = dom["dune_rate_m_yr"].to_numpy(dtype=float)

    y_tick = 2.0
    half = float(np.nanmax(np.abs(np.concatenate([lrr, ept, dune])))) + 1.0
    half = float(np.ceil(half / y_tick) * y_tick)
    ax.set_xlim(GIS_FIRST - 0.5, GIS_LAST + 0.5)
    ax.set_ylim(-half, half)
    town_bands(ax)

    ax.axhline(0, color=INK_MUTED, lw=0.6, zorder=2)
    ax.plot(g, ept, color=C_ENDPOINT, lw=0.9, ls=(0, (3, 1.5)), zorder=3)
    ax.plot(g, lrr, color=C_LRR, lw=1.0, zorder=4)
    ax.plot(g, dune, color=C_DUNE, lw=1.1, zorder=5)

    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_minor_locator(MultipleLocator(5))
    ax.yaxis.set_major_locator(MultipleLocator(y_tick))
    ax.yaxis.grid(True, zorder=0)
    ax.set_axisbelow(True)
    open_frame(ax)
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    ax.set_ylabel("Rate of change (m/yr)")

    # The legend sits outside the axes and shrinks them at draw time, so it
    # goes on BEFORE the structure labels are placed against the data.
    handles = [
        Line2D([], [], color=C_DUNE, lw=1.1, label="Dune line (two surveys)"),
        Line2D([], [], color=C_LRR, lw=1.0,
               label="CoastSat shoreline, LRR"),
        Line2D([], [], color=C_ENDPOINT, lw=0.9, ls=(0, (3, 1.5)),
               label="CoastSat shoreline, endpoint"),
    ]
    fig.legend(handles=handles, loc="outside upper right", ncol=3)
    structures(ax)

    caption(fig, (
        f"Alongshore rate of change {start}-{end} by GIS domain, seaward "
        f"positive. Red: the digitized dune line, its {start} and {end} "
        f"positions differenced over the survey interval. Dark blue: the "
        f"CoastSat shoreline as the per-domain linear regression rate over "
        f"the window; light blue, dashed: the CoastSat endpoint rate from the "
        f"mean position in a one-year window about each survey date. Where "
        f"the two blues part is where the choice of estimator matters. Red and "
        f"blue here mark the two features, not the two vintages. "
        f"Bands mark Buxton, Avon and the Tri-Village; "
        f"the solid hairline is the Buxton groin, the dotted ones the Avon "
        f"and Rodanthe piers. Domain 1 is Cape Point, 90 is Pea Island."))
    save(fig, out / "alongshore_dune_vs_coastsat", close=True)


# -----------------------------------------------------------------------------
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="dune-line change against the CoastSat shoreline")
    ap.add_argument("--start-year", type=int, default=1984)
    ap.add_argument("--end-year", type=int, default=2004)
    ap.add_argument("--start-date", default=None,
                    help="survey date of the start line, YYYY-MM-DD")
    ap.add_argument("--end-date", default=None,
                    help="survey date of the end line, YYYY-MM-DD")
    ap.add_argument("--half-window-days", type=float, default=SIX_MONTHS_DAYS,
                    help="half-width of the CoastSat window about each survey "
                         "date (default six months)")
    a = ap.parse_args(argv)

    assumed = {}
    dates = {}
    for key, yr in (("start", a.start_year), ("end", a.end_year)):
        given = getattr(a, f"{key}_date") or KNOWN_SURVEY_DATES.get(yr)
        if given is None:
            given = f"{yr}-07-01"
            assumed[key] = True
        dates[key] = datetime.strptime(given, "%Y-%m-%d").replace(tzinfo=timezone.utc)
    d0, d1 = dates["start"], dates["end"]
    years = (d1 - d0).days / DAYS_PER_YEAR
    print(f"survey dates  {d0.date()} -> {d1.date()}  ({years:.2f} yr)"
          + ("   [END DATE ASSUMED]" if "end" in assumed else "")
          + ("   [START DATE ASSUMED]" if "start" in assumed else ""))

    apply_style()
    out = OUT_ROOT / f"{a.start_year}_{a.end_year}"
    out.mkdir(parents=True, exist_ok=True)
    sup = support_dir(out)      # figures at the top, everything else here

    # dune
    p0 = dune_position_by_domain(a.start_year)
    p1 = dune_position_by_domain(a.end_year)
    dune_change = -(p1 - p0)               # seaward positive
    dune_rate = dune_change / years
    print(f"dune domains with both years: {int((p0.notna() & p1.notna()).sum())} / 90")

    # coastsat
    lookup = pd.read_csv(transect_lookup())
    lookup = lookup[lookup["domain_number"].between(GIS_FIRST, GIS_LAST)]
    cache: dict = {}
    ep = endpoint_by_transect(lookup, d0, d1, a.half_window_days, cache)
    for c in ("first_obs_start", "last_obs_start", "first_obs_end", "last_obs_end"):
        ep[c] = pd.to_datetime(ep[c], utc=True).dt.strftime("%Y-%m-%d")
    ep.to_csv(sup / "transect_coastsat_endpoint.csv", index=False,
              float_format="%.3f")
    # a window truncated by the record itself is not a full seasonal cycle
    span = {}
    for key, d in (("start", d0), ("end", d1)):
        f = pd.to_datetime(ep[f"first_obs_{key}"]).median()
        l = pd.to_datetime(ep[f"last_obs_{key}"]).median()
        lo = (d - timedelta(days=a.half_window_days)).date()
        hi = (d + timedelta(days=a.half_window_days)).date()
        span[key] = (lo, hi, f.date(), l.date(),
                     (f.date() - lo).days > 45 or (hi - l.date()).days > 45)
    lrr = pd.read_csv(lrr_csv(a.start_year, a.end_year))
    lrr = lrr[lrr["domain_number"].between(GIS_FIRST, GIS_LAST)]

    ep_dom = ep.groupby("domain_number").agg(
        cs_endpoint_m_yr=("endpoint_rate_m_yr", "mean"),
        cs_endpoint_change_m=("change_m", "mean"),
        n_obs_start=("n_start", "median"),
        n_obs_end=("n_end", "median"),
        n_transects=("transect_id", "size"))
    lrr_dom = lrr.groupby("domain_number").agg(cs_lrr_m_yr=("lrr_m_yr", "mean"))

    dom = pd.DataFrame({"gis": range(GIS_FIRST, GIS_LAST + 1)}).set_index("gis")
    dom["dune_start_m"] = p0
    dom["dune_end_m"] = p1
    dom["dune_change_m"] = dune_change
    dom["dune_rate_m_yr"] = dune_rate
    dom = dom.join(lrr_dom.rename_axis("gis")).join(ep_dom.rename_axis("gis"))
    dom["dune_minus_lrr_m_yr"] = dom["dune_rate_m_yr"] - dom["cs_lrr_m_yr"]
    dom["dune_minus_endpoint_m_yr"] = dom["dune_rate_m_yr"] - dom["cs_endpoint_m_yr"]
    dom = dom.reset_index()
    dom.to_csv(sup / "domain_comparison.csv", index=False, float_format="%.3f")

    stats = {c: _fit_stats(dom["dune_rate_m_yr"].to_numpy(), dom[c].to_numpy())
             for c in ("cs_lrr_m_yr", "cs_endpoint_m_yr")}
    both = _fit_stats(dom["cs_lrr_m_yr"].to_numpy(), dom["cs_endpoint_m_yr"].to_numpy())

    # sensitivity of the endpoint rate to the window centre, for an assumed date
    sens = []
    if assumed:
        for label, shift in (("-6 mo", -SIX_MONTHS_DAYS), ("0", 0.0),
                             ("+6 mo", SIX_MONTHS_DAYS)):
            dd0 = d0 + timedelta(days=shift if "start" in assumed else 0)
            dd1 = d1 + timedelta(days=shift if "end" in assumed else 0)
            e = endpoint_by_transect(lookup, dd0, dd1, a.half_window_days, cache)
            e_dom = e.groupby("domain_number")["endpoint_rate_m_yr"].mean()
            st = _fit_stats(
                dom.set_index("gis")["dune_rate_m_yr"].reindex(e_dom.index).to_numpy(),
                e_dom.to_numpy())
            sens.append((label, dd0.date(), dd1.date(), float(e_dom.mean()), st))

    scatter_figure(dom, out, a.start_year, a.end_year, stats)
    alongshore_figure(dom, out, a.start_year, a.end_year)

    # provenance
    isl = dom[["dune_rate_m_yr", "cs_lrr_m_yr", "cs_endpoint_m_yr"]].mean()
    known_src = {
        1984: "USGS 1984 aerial photo, Henderson release "
              "(`D:\\Hatteras_GIS\\Aerial\\1984_henderson\\1984_metadata`)",
        2004: "Google Earth capture date (the raw_GE frames carry none); "
              "Hannah, 2026-09-15",
    }

    def _src(key, yr):
        if key in assumed:
            return "**ASSUMED mid-year**; no date known for this line"
        if getattr(a, f"{key}_date"):
            return "given on the command line"
        return known_src.get(yr, "KNOWN_SURVEY_DATES in the script")

    start_src = _src("start", a.start_year)
    end_src = _src("end", a.end_year)
    lines = [
        f"# Dune line vs CoastSat shoreline, {a.start_year}-{a.end_year}",
        "",
        f"Written {datetime.now():%Y-%m-%d %H:%M} by "
        f"`scripts/input_prep/5-scr/duneline_vs_coastsat/duneline_vs_coastsat.py`.",
        "",
        "## Inputs",
        "",
        f"* dune lines: `2-brie-offset/raw_offsets/{a.start_year}_duneline_offset_raw.csv`, "
        f"`{a.end_year}_duneline_offset_raw.csv` (first row per transect, domain mean, "
        f"as `hindcast.load_absolute_dune_distance`). Both are ArcGIS exports, so the "
        f"1 m GIS-vs-shapely convention cancels.",
        f"* CoastSat LRR: `{lrr_csv(a.start_year, a.end_year).relative_to(SCR_ROOT).as_posix()}` "
        f"(window {a.start_year}-01-01 to {a.end_year}-12-31, per-transect OLS).",
        f"* CoastSat endpoint: mean chainage within ±{a.half_window_days:.0f} days of each "
        f"survey date, per transect, from `coastsat_timeseries/`.",
        "* transect → domain: `transect_domains/transect_domain_lookup.csv`.",
        "",
        "## Survey dates",
        "",
        "| line | date | source |",
        "|---|---|---|",
        f"| {a.start_year} | {d0.date()} | {start_src} |",
        f"| {a.end_year} | {d1.date()} | {end_src} |",
        "",
        f"Survey interval {years:.2f} yr. Sign: seaward positive in every column; "
        f"a negative rate is retreat. Dune change is `-(ORIG_LEN_end - ORIG_LEN_start)`.",
        "",
        "## Island-wide means (m/yr)",
        "",
        "| dune line | CoastSat LRR | CoastSat endpoint |",
        "|---|---|---|",
        f"| {isl['dune_rate_m_yr']:.2f} | {isl['cs_lrr_m_yr']:.2f} | {isl['cs_endpoint_m_yr']:.2f} |",
        "",
        "## Agreement, per domain (y against dune rate x)",
        "",
        "| y | n | r | slope | intercept | RMSE | bias (y − x) |",
        "|---|---|---|---|---|---|---|",
    ]
    for c, name in (("cs_lrr_m_yr", "CoastSat LRR"), ("cs_endpoint_m_yr", "CoastSat endpoint")):
        s = stats[c]
        lines.append(f"| {name} | {s['n']} | {s['r']:.2f} | {s['slope']:.2f} | "
                     f"{s['intercept']:.2f} | {s['rmse']:.2f} | {s['bias']:+.2f} |")
    lines += [
        "",
        f"CoastSat endpoint against CoastSat LRR (two estimators of the same "
        f"series): r = {both['r']:.2f}, slope = {both['slope']:.2f}, RMSE = {both['rmse']:.2f}, "
        f"bias = {both['bias']:+.2f} m/yr.",
        "",
        "## Window occupancy",
        "",
        f"Median CoastSat observations per transect inside the start window: "
        f"{int(ep['n_start'].median())}; inside the end window: {int(ep['n_end'].median())}. "
        f"Transects with an empty window: start {int((ep['n_start'] == 0).sum())}, "
        f"end {int((ep['n_end'] == 0).sum())}, of {len(ep)}.",
        "",
        "| window | asked for | median first obs | median last obs | truncated |",
        "|---|---|---|---|---|",
    ]
    for key in ("start", "end"):
        lo, hi, f, l, trunc = span[key]
        lines.append(f"| {key} | {lo} to {hi} | {f} | {l} | "
                     f"{'**yes**' if trunc else 'no'} |")
    if any(v[4] for v in span.values()):
        lines += [
            "",
            "A truncated window does not average a full seasonal cycle. The "
            "CoastSat record begins 1984-09-21 on most transects (1984-05/06 on "
            "a few), so a window centred on the 1984-09-19 photo holds only the "
            "autumn and winter after it.",
        ]
    if sens:
        lines += [
            "",
            "## Sensitivity of the endpoint rate to the assumed survey date",
            "",
            "| centre shift | start | end | island mean endpoint (m/yr) | r vs dune | slope | RMSE |",
            "|---|---|---|---|---|---|---|",
        ]
        for label, s0, s1, m, st in sens:
            lines.append(f"| {label} | {s0} | {s1} | {m:.2f} | {st['r']:.2f} | "
                         f"{st['slope']:.2f} | {st['rmse']:.2f} |")
    lines += [
        "",
        "## Read this before quoting it",
        "",
        "* The dune line and the waterline are different features. A gap between "
        "them is beach-width change as much as it is disagreement.",
        "* The LRR spans the calendar window; the endpoint spans the survey interval. "
        "They are not the same length of record.",
        "* `n_obs_*` is the median per-transect count inside a one-year window. "
        "One storm inside a window moves that end.",
    ]
    (sup / "PROVENANCE.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"island means  dune {isl['dune_rate_m_yr']:+.2f}  lrr {isl['cs_lrr_m_yr']:+.2f}  "
          f"endpoint {isl['cs_endpoint_m_yr']:+.2f} m/yr")
    for c in stats:
        s = stats[c]
        print(f"  {c:18s} n={s['n']} r={s['r']:.2f} slope={s['slope']:.2f} "
              f"rmse={s['rmse']:.2f} bias={s['bias']:+.2f}")
    print(f"  endpoint vs lrr    r={both['r']:.2f} slope={both['slope']:.2f} bias={both['bias']:+.2f}")
    for label, s0, s1, m, st in sens:
        print(f"  sens {label:6s} {s0} {s1}  mean={m:+.2f}  r={st['r']:.2f}")
    print(f"-> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
