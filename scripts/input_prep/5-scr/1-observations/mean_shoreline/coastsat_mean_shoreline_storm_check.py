"""
Was a window mean shaped by a storm? A check on the +/-1 yr means that become the shoreline offset.

    python scripts/input_prep/5-scr/1-observations/mean_shoreline/coastsat_mean_shoreline_storm_check.py
    python scripts/input_prep/5-scr/1-observations/mean_shoreline/coastsat_mean_shoreline_storm_check.py --centred-on alace_1996

Storms in and around the window, the island-wide anomaly through it, and
the shift post-storm passes make; writes a figure and README per window. Details: scripts/input_prep/5-scr/1-observations/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-29
"""

from __future__ import annotations

import argparse
import sys
from datetime import date
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_REPO = next(_p for _p in Path(__file__).resolve().parents
             if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(_REPO / "scripts"))
sys.path.insert(0, str(_REPO / "scripts" / "input_prep" / "5-scr" / "lib"))
sys.path.insert(0, str(_REPO / "scripts" / "input_prep" / "3-env-forcings" / "3-storms"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import scr_paths  # noqa: E402,F401

import storm_figures as sf  # noqa: E402
from coastsat_lrr import load_timeseries  # noqa: E402
from coastsat_mean_shoreline import (  # noqa: E402
    SURVEY_ANCHORS, Window, dune_line_date, timeseries_file,
)
from site_layer import hat_env_forcings as env  # noqa: E402
from site_layer import hat_figure_style as fs  # noqa: E402
from site_layer.hat_observed_rates import mean_shoreline_csv, mean_shoreline_dir  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
RECORD = (1984, 2024)
DEFAULT_CONTEXT_YEARS = 3
# Beach recovery after a storm runs weeks to months
DEFAULT_RECOVERY_DAYS = 90
DEFAULT_MIN_COVERAGE = 0.5
MIN_OBS = 10          # the producer's minimum for a transect mean
CHECK_DIR = "storm_check"
# -----------------------------------------------------------------------------


# The storm record, 1984-2024

# Every event of 1984-2024, one row each, with its peak hour, Rhigh (m MHW), hours above the berm ...
def load_events(variant):
    forcing = sf.load_forcing()
    parts = []
    for (a, b), keep in (((1984, 2004), lambda y: y <= 2003), ((2004, 2024), lambda y: y >= 2004)):
        df = pd.read_csv(env.storm_summary_file(a, b, variant), parse_dates=["StartTime", "EndTime"])
        parts.append(df[keep(df.calendar_year)])
    df = pd.concat(parts, ignore_index=True).rename(columns={"period": "tp_s"})
    df["peak"] = [forcing.TWL.loc[s:e].idxmax() for s, e in zip(df.StartTime, df.EndTime)]
    df["rhigh_m"] = df.Rhigh * 10.0
    rebuilt = forcing.TWL.loc[df.peak].values - sf.MHW_NAVD88
    if not np.allclose(rebuilt, df.rhigh_m, atol=0.01):
        raise RuntimeError("rebuilt TWL does not reproduce the summary's Rhigh; the peak hours are wrong")
    trimmed = df["trimmed_from"] if "trimmed_from" in df else pd.Series(0, index=df.index)
    df["raw_hours"] = np.where(trimmed > 0, trimmed, df.duration).astype(int)

    # storm type, as storm_figures.classify, but from 1983 so 1984 is covered
    h = sf.load_hurdat2(first_year=1983)
    win = pd.Timedelta(hours=sf.TC_WINDOW_H)
    names, kms = [], []
    for t in df.peak:
        c = h[(h.t >= t - win) & (h.t <= t + win)]
        r = c.loc[c.km.idxmin()] if not c.empty else None
        names.append(r.tc_name if r is not None else None)
        kms.append(r.km if r is not None else np.nan)
    df["tc_name"], df["tc_km"] = names, kms
    df["type"] = np.where(df.tc_km <= sf.TC_NEAR_KM, "tropical", "other")

    df["rank_rhigh"] = df.rhigh_m.rank(ascending=False, method="min").astype(int)
    df["pct_rhigh"] = 100.0 * df.rhigh_m.rank(pct=True)
    df["rank_hours"] = df.raw_hours.rank(ascending=False, method="min").astype(int)
    return df.sort_values("peak").reset_index(drop=True)


# Median over the record of each year's highest event
def annual_max_median(events, col, empty):
    ann = events.groupby(events.peak.dt.year)[col].max()
    ann = ann.reindex(range(RECORD[0], RECORD[1] + 1)).fillna(empty)   # a year with no event
    return float(ann.median())


# 'height', 'length', 'both' or '' per event
def major_by(df, major):
    h, l = df.rhigh_m >= major["rhigh"], df.raw_hours >= major["hours"]
    return np.select([h & l, h, l], ["both", "height", "length"], "")


# A tropical event's name, or None
def event_name(r):
    if r["type"] == "tropical":
        return sf.label_text(pd.Series(dict(name=r.tc_name, peak=r.peak, tc_status="")))
    return None


# Storm-hours above the berm in every 2-yr span of the record, stepped monthly
def storm_hours_2yr(events, step_days=30):
    lo, hi = pd.Timestamp(f"{RECORD[0]}-01-01"), pd.Timestamp(f"{RECORD[1]}-12-31")
    rows = []
    t = lo
    while t + pd.DateOffset(years=2) <= hi + pd.Timedelta(days=1):
        e = t + pd.DateOffset(years=2)
        sel = events[(events.peak >= t) & (events.peak < e)]
        rows.append((t, e, int(sel.raw_hours.sum()), len(sel),
                     float(sel.rhigh_m.max()) if len(sel) else np.nan))
        t = t + pd.Timedelta(days=step_days)
    return pd.DataFrame(rows, columns=["start", "end", "hours", "events", "max_rhigh"])


# The shoreline

# Every CoastSat position of the window's included transects over the context span, as the anomaly ...
def load_positions(window, context_lo, context_hi):
    means = pd.read_csv(mean_shoreline_csv(*window.key))
    means = means[means.included]
    parts = []
    for tid, dom, m in zip(means.transect_id, means.domain_number, means.mean_chainage_m):
        ts = load_timeseries(str(timeseries_file(tid)))
        ts["date"] = ts["date"].dt.tz_localize(None)
        ts = ts[(ts.date >= context_lo) & (ts.date <= context_hi + pd.Timedelta(days=1))]
        parts.append(ts.assign(transect_id=tid, domain_number=int(dom),
                               window_mean=m, anomaly=ts.chainage_m - m))
    return pd.concat(parts, ignore_index=True), means


# Per-day island median anomaly and transect count
def island_series(pos, n_transects, min_coverage):
    day = pos.assign(day=pos.date.dt.normalize()).groupby("day")
    s = pd.DataFrame({"n_transects": day.transect_id.nunique(),
                      "median_anomaly_m": day.anomaly.median(),
                      "q25_m": day.anomaly.quantile(0.25),
                      "q75_m": day.anomaly.quantile(0.75)})
    s["coverage"] = s.n_transects / n_transects
    s["drawn"] = s.coverage >= min_coverage
    return s.reset_index().rename(columns={"day": "date"})


# Each transect's window mean without the passes inside `recovery_days` after any of `storms`
def mean_shift(pos, window, storms, recovery_days):
    inw = pos[(pos.date >= window.lo) & (pos.date <= window.hi + pd.Timedelta(days=1))]
    rec = pd.Timedelta(days=recovery_days)
    g = inw.groupby("transect_id")
    out = pd.DataFrame({"domain_number": g.domain_number.first(),
                        "n_obs": g.size(), "window_mean_m": g.chainage_m.mean()})
    masks = {}
    for _, s in storms.iterrows():
        masks[s.key] = (inw.date >= s.peak) & (inw.date < s.peak + rec)
    if len(masks) > 1:
        masks["all_major"] = np.logical_or.reduce(list(masks.values()))
    for key, m in masks.items():
        kept = inw[~m].groupby("transect_id").chainage_m
        out[f"n_after_{key}"] = inw[m].groupby("transect_id").size().reindex(out.index, fill_value=0)
        out[f"mean_without_{key}_m"] = kept.mean().reindex(out.index)
        out[f"shift_{key}_m"] = out[f"mean_without_{key}_m"] - out.window_mean_m
        short = kept.size().reindex(out.index, fill_value=0) < MIN_OBS
        out.loc[short, f"shift_{key}_m"] = np.nan
    return out.reset_index(), list(masks)


# Figure

# Shade the window
def _shade(ax, window):
    ax.axvspan(window.lo, window.hi, color=fs.C["BASE_FILL"], alpha=0.55, lw=0, zorder=0)


# The storms around the window, one panel
def figure(ctx, spans, window, major, folder, label, anchor, dune):
    fig, a = plt.subplots(figsize=fs.figsize("double", height=3.4))
    xlim = (pd.Timestamp(window.lo) - pd.DateOffset(years=spans["ctx"]),
            pd.Timestamp(window.hi) + pd.DateOffset(years=spans["ctx"]))
    halo = fs._halo(2.2)

    _shade(a, window)
    if anchor:
        f0, f1 = (pd.Timestamp(d) for d in SURVEY_ANCHORS[anchor]["flown"])
        a.axvline(f0 + (f1 - f0) / 2, color=fs.INK, lw=0.9, zorder=2)
    if dune:
        a.axvline(pd.Timestamp(dune[1]), color=fs.C_1984, lw=0.9, ls=(0, (1, 1.5)), zorder=2)
    a.text(window.lo + (window.hi - window.lo) / 2, 0.985, "averaging window",
           transform=a.get_xaxis_transform(), ha="center", va="top", fontsize=7,
           color=fs.INK_MUTED, path_effects=halo, zorder=7)

    for k, size in (("other", 14), ("tropical", 20)):
        d = ctx[ctx.type == k]
        a.scatter(d.peak, d.rhigh_m, s=size, color=sf.TYPE_COLOURS[k], edgecolors="white",
                  linewidths=0.35, zorder=3)
    big = ctx[ctx.major]
    a.scatter(big.peak, big.rhigh_m, s=46, facecolors="none", edgecolors=fs.INK, linewidths=0.8, zorder=4)
    a.axhline(major["rhigh"], color=fs.INK, lw=0.8, ls=(0, (5, 2)), zorder=1)
    a.text(xlim[1], major["rhigh"] - 0.03, f"median annual maximum 1984–2024 ({major['rhigh']:.2f} m)  ",
           ha="right", va="top", fontsize=6.8, color=fs.INK, path_effects=halo)
    # Named storms carry their HURDAT2 name
    for _, r in big.iterrows():
        txt = r["name"] or (f"{r.raw_hours} h above berm" if r.major_by == "length" else None)
        if txt:
            a.annotate(txt, (r.peak, r.rhigh_m), xytext=(0, 5), textcoords="offset points",
                       ha="center", va="bottom", fontsize=6.5, color=fs.INK, path_effects=halo, zorder=6)
    a.set_xlim(*xlim)
    a.set_ylabel(r"$R_{high}$ (m above MHW)")
    a.set_ylim(sf.BERM_MHW - 0.05, max(ctx.rhigh_m.max(), major["rhigh"]) + 0.35)
    a.yaxis.grid(True, zorder=0)
    fs.open_frame(a)

    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    handles = [Patch(color=sf.TYPE_COLOURS[k], label=sf.TYPE_LABELS[k]) for k in sf.TYPES]
    handles += [Line2D([], [], marker="o", ls="", mfc="none", mec=fs.INK, label="major storm")]
    if anchor:
        handles.append(Line2D([], [], color=fs.INK, lw=0.9, label="lidar flights (window centre)"))
    if dune:
        handles.append(Line2D([], [], color=fs.C_1984, lw=0.9, ls=(0, (1, 1.5)),
                              label=f"{dune[0]} dune-line imagery"))
    a.legend(handles=handles, loc="lower center", ncol=3, bbox_to_anchor=(0.5, 1.02),
             fontsize=7, handlelength=1.4, frameon=False)

    fs.caption(fig, (
        f"Storms around the {window.span} CoastSat mean shoreline ({window.described}). "
        f"Grey band: the averaging window. Every high-water event of the hindcast storm series "
        f"({spans['variant']}) from {spans['ctx']} yr before the window to {spans['ctx']} yr after, at its "
        f"peak total water level (Duck gauge + Stockdon 2006 R2% from WIS 63228 waves; an estimate, not "
        f"an observation at Hatteras). Circled: major, Rhigh or hours above the berm at or above its "
        f"median annual maximum of 1984–2024 ({major['rhigh']:.2f} m MHW, dashed; {major['hours']:.0f} h), "
        f"a level reached in half the years. Tropical cyclones within {sf.TC_NEAR_KM} km are labelled with "
        f"their HURDAT2 name; an unnamed storm major by length only is labelled with its hours above the "
        f"berm. Each storm's rank in the {spans['n_record']} events of 1984–2024 is in the README and "
        f"events table beside this figure."))
    fs.save(fig, folder / f"storm_check_{label}.png", bbox_inches="tight")
    plt.close(fig)


# The readme

# A signed number, or n/a
def _pct(v):
    return f"{v:+.1f}" if np.isfinite(v) else "n/a"


# The README beside one window's check
def write_readme(folder, label, window, anchor, dune, ctx, major, n_record, spans_df,
                 window_hours, shift, keys, storms, recovery_days, series, n_transects, variant, ctx_years):
    inw = ctx[(ctx.peak >= window.lo) & (ctx.peak <= window.hi + pd.Timedelta(days=1))]
    pct_hours = 100.0 * (spans_df.hours < window_hours).mean()
    rank_hours = int((spans_df.hours > window_hours).sum()) + 1

    def ev_row(r):
        where = ("before" if r.peak < window.lo else "after" if r.peak > window.hi else "inside")
        return (f"| {r.peak:%Y-%m-%d} | {r['name'] or '—'} | {r.type} | {r.rhigh_m:.2f} | "
                f"{r.rank_rhigh} | {r.raw_hours} | {r.rank_hours} | {r.major_by or '—'} | {where} |")

    big = ctx[ctx.major]
    top_hours = ctx.nlargest(5, "raw_hours")
    big_rows = "\n".join(ev_row(r) for _, r in big.iterrows()) or "| none | | | | | | | | |"
    hour_rows = "\n".join(ev_row(r) for _, r in top_hours.sort_values("peak").iterrows())

    shift_rows = []
    for key in keys:
        col = shift[f"shift_{key}_m"]
        n_aft = shift[f"n_after_{key}"]
        nm = "all major storms together" if key == "all_major" else key
        shift_rows.append(f"| {nm} | {int(n_aft.median())} | {_pct(col.median())} | "
                          f"{_pct(col.quantile(0.05))} to {_pct(col.quantile(0.95))} | "
                          f"{_pct(col.abs().max())} | {int(col.isna().sum())} |")
    shift_block = ("| storm (peak) | passes removed per transect (median) | shift, median (m) | "
                   "shift, 5th to 95th pct (m) | largest |shift| (m) | transects left < 10 passes |\n"
                   "|---|---|---|---|---|---|\n" + "\n".join(shift_rows)) if keys else \
        "No major storm peaked inside the window or in the recovery period before it, so nothing was removed."

    drawn = series[series.drawn]
    inwin = drawn[(drawn.date >= window.lo) & (drawn.date <= window.hi)]
    # the one-paragraph answer, computed
    if len(storms):
        allk = "all_major" if "all_major" in keys else keys[0]
        med = shift[f"shift_{allk}_m"].median()
        verdict = (f"{len(storms)} major storm(s) peaked inside the window or within {recovery_days} days "
                   f"before it. Dropping the passes in the {recovery_days} days after "
                   f"{'them' if len(storms) > 1 else 'it'} moves the island-median transect mean by "
                   f"**{med:+.1f} m** (positive: the mean would sit further seaward without the post-storm "
                   f"passes), against a median standard error of a transect mean of "
                   f"~{spans_df.attrs.get('se_med', float('nan')):.1f} m and the 10 m Barrier3D cell.")
    else:
        verdict = "No major storm peaked inside the window or within the recovery period before it."
    verdict += (f" The window held **{window_hours} storm-hours** above the berm; of the 2-yr spans of "
                f"{RECORD[0]}–{RECORD[1]} it ranks {rank_hours} of {len(spans_df)} "
                f"(stormier than {pct_hours:.0f}% of them).")

    text = f"""# storm_check/{label} -- were there big storms around this window mean?

Written by `scripts/input_prep/5-scr/1-observations/mean_shoreline/coastsat_mean_shoreline_storm_check.py`
on {date.today().isoformat()}. A check, not an input: nothing here changes
the mean line in the folder above.

The mean line is {window.described}. The question is whether a
big storm inside that window, or just before it, pulled the mean landward.

## Answer

{verdict}

## 1. The major storms, {ctx_years} yr either side of the window

**Major** = at or above the median annual maximum of {RECORD[0]}–{RECORD[1]} in
**height** (Rhigh ≥ {major['rhigh']:.2f} m above MHW, {major['rhigh'] + sf.MHW_NAVD88:.2f} m NAVD88) or in
**length** (≥ {major['hours']:.0f} h above the berm): a level the record reaches in half its
years. Ranks are among the {n_record} events of {RECORD[0]}–{RECORD[1]} (1 = highest / longest).
Hours are above the berm before the 24 h trim.

| peak | name | type | Rhigh (m MHW) | rank, Rhigh | hours | rank, hours | major by | vs window |
|---|---|---|---|---|---|---|---|---|
{big_rows}

The five longest events in the same span (a long nor'easter can move more
sand than a higher, shorter storm):

| peak | name | type | Rhigh (m MHW) | rank, Rhigh | hours | rank, hours | major by | vs window |
|---|---|---|---|---|---|---|---|---|
{hour_rows}

Inside the window itself: {len(inw)} events, {int(inw.major.sum())} of them major;
highest {inw.rhigh_m.max():.2f} m ({'' if inw.empty else inw.loc[inw.rhigh_m.idxmax(), 'peak'].strftime('%Y-%m-%d')}).

## 2. Was the window stormy?

{window_hours} storm-hours above the berm fell inside the window. Over every
2-yr span of {RECORD[0]}–{RECORD[1]} (stepped monthly) the median is
{int(spans_df.hours.median())} h and the range {int(spans_df.hours.min())}–{int(spans_df.hours.max())} h,
so this window ranks {rank_hours} of {len(spans_df)}.

## 3. Did a storm move the mean?

Each transect's window mean recomputed without the CoastSat passes that
fall within **{recovery_days} days after** a major storm's peak (storms
peaking inside the window, or within {recovery_days} days before it). The
shift is (mean without) − (mean as built); positive means the post-storm
passes had pulled the mean **landward**. A transect left with fewer than
{MIN_OBS} passes gets no shift.

{shift_block}

Per transect in `storm_check_{label}_mean_shift.csv`.

## 4. The shoreline through the window

The island-median CoastSat position, each transect relative to its
own window mean, from {ctx_years} yr before to {ctx_years} yr after. {len(drawn)} image
dates have at least {int(100 * series.attrs['min_cov'])}% of the {n_transects} transects
({len(inwin)} of them inside the window). A storm that moved the shoreline
shows as a step down after its peak that does not recover.
Series in `supporting/storm_check_{label}_island_series.csv`.

## Caveats

- The storm water levels are the model's estimate (Duck gauge + Stockdon
  R2% from WIS waves, beach slope 0.06), not observations at Hatteras.
- "Major" and the {recovery_days}-day recovery are choices; re-run with
  `--major-rhigh`, `--major-hours` or `--recovery-days` to test them.
- The shift is the storm's effect *through the passes that followed it*. A
  storm that moved the shoreline for good moves every later pass too, and
  that part is not removable by dropping passes: look at the island series for it.
- Landsat 5 alone before 1999, so the 1996 window has fewer passes to drop.

## Files

| file | what it is |
|---|---|
| `storm_check_{label}.png` | every storm around the window, major ones circled |
| `storm_check_{label}_events.csv` | every event in the context span: peak, Rhigh, hours, type, HURDAT2 name and distance, ranks in 1984–2024, major, inside window |
| `storm_check_{label}_mean_shift.csv` | per transect: window mean, passes removed, mean without, shift |
| `supporting/storm_check_{label}_island_series.csv` | per image date: transects, median and quartiles of the anomaly, coverage |
| `supporting/storm_check_{label}.pdf`, `supporting/CAPTIONS.md` | vector figure, caption |

Storm series: `{variant}` (`hat_env_forcings.DEFAULT_STORM_VARIANT`).
"""
    (folder / "README.md").write_text(text, encoding="utf-8")


# One row in the window's PROVENANCE.md Files table, if it is missing
def link_from_provenance(window_dir):
    p = window_dir / "PROVENANCE.md"
    if not p.is_file():
        return
    text = p.read_text(encoding="utf-8")
    if "`storm_check/`" in text:
        return
    row = ("| `storm_check/` | were there big storms around this window? The storm record "
           "3 yr either side, ranked in 1984-2024, and the mean without post-storm passes; "
           "written by `coastsat_mean_shoreline_storm_check.py` |\n")
    tail = "\n\nResolved through `hat_observed_rates"
    if tail in text:
        p.write_text(text.replace(tail, "\n" + row + tail[1:], 1), encoding="utf-8")


# One window: load, compute, draw, write
def run(window, events, major, variant, ctx_years, recovery_days, min_cov):
    label = window.label
    window_dir = mean_shoreline_dir(*window.key)
    if not mean_shoreline_csv(*window.key).is_file():
        raise FileNotFoundError(f"no mean shoreline for {label}; run coastsat_mean_shoreline.py first")
    folder = window_dir / CHECK_DIR
    folder.mkdir(exist_ok=True)

    lo = window.lo - pd.DateOffset(years=ctx_years)
    hi = window.hi + pd.DateOffset(years=ctx_years)
    ctx = events[(events.peak >= lo) & (events.peak <= hi + pd.Timedelta(days=1))].copy()
    ctx["name"] = [event_name(r) for _, r in ctx.iterrows()]
    ctx["major_by"] = major_by(ctx, major)
    ctx["major"] = ctx.major_by != ''
    ctx["vs_window"] = np.where(ctx.peak < window.lo, "before",
                                np.where(ctx.peak > window.hi + pd.Timedelta(days=1), "after", "inside"))

    pos, means = load_positions(window, lo, hi)
    n_tr = len(means)
    series = island_series(pos, n_tr, min_cov)
    series.attrs["min_cov"] = min_cov

    rec = pd.Timedelta(days=recovery_days)
    storms = ctx[ctx.major & (ctx.peak >= window.lo - rec) & (ctx.peak <= window.hi)].copy()
    storms["key"] = [f"{r.name or 'event'}_{r.peak:%Y-%m-%d}".replace(" ", "_")
                     for r in storms.itertuples()]
    shift, keys = mean_shift(pos, window, storms, recovery_days)

    spans_df = storm_hours_2yr(events)
    spans_df.attrs["se_med"] = float(means.se_chainage_m.median())
    wsel = events[(events.peak >= window.lo) & (events.peak <= window.hi + pd.Timedelta(days=1))]
    window_hours = int(wsel.raw_hours.sum())

    dune = dune_line_date(window.period_start)
    figure(ctx, dict(ctx=ctx_years, variant=variant, n_record=len(events)),
           window, major, folder, label, window.anchor, dune)

    cols = ["peak", "rhigh_m", "raw_hours", "tp_s", "type", "name", "tc_name", "tc_km",
            "rank_rhigh", "pct_rhigh", "rank_hours", "major", "major_by", "vs_window"]
    out = ctx[cols].rename(columns={"rhigh_m": "rhigh_m_mhw", "raw_hours": "hours_above_berm",
                                    "tc_name": "nearest_tc", "tc_km": "nearest_tc_km"})
    out["nearest_tc_km"] = out.nearest_tc_km.round(0)
    out.round(3).to_csv(folder / f"storm_check_{label}_events.csv", index=False)
    shift.round(3).to_csv(folder / f"storm_check_{label}_mean_shift.csv", index=False)
    series.round(3).to_csv(fs.support_dir(folder) / f"storm_check_{label}_island_series.csv", index=False)
    write_readme(folder, label, window, window.anchor, dune, ctx, major, len(events), spans_df,
                 window_hours, shift, keys, storms, recovery_days, series, n_tr, variant, ctx_years)
    link_from_provenance(window_dir)

    print(f"\n== {label} -> {folder}")
    print(ctx[ctx.major][["peak", "name", "type", "rhigh_m", "rank_rhigh", "raw_hours", "major_by", "vs_window"]]
          .to_string(index=False))
    print(f"window storm-hours {window_hours}; 2-yr median {int(spans_df.hours.median())}")
    for k in keys:
        c = shift[f"shift_{k}_m"]
        print(f"  shift without {k}: median {c.median():+.2f} m, p5 {c.quantile(.05):+.2f}, p95 {c.quantile(.95):+.2f}")


# Run: the chosen windows
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--centred-on", choices=sorted(SURVEY_ANCHORS), action="append")
    g.add_argument("--window-dates", nargs=2, metavar=("FIRST", "LAST"))
    g.add_argument("--window", nargs=2, type=int, metavar=("START", "END"))
    ap.add_argument("--context-years", type=int, default=DEFAULT_CONTEXT_YEARS)
    ap.add_argument("--recovery-days", type=int, default=DEFAULT_RECOVERY_DAYS)
    ap.add_argument("--min-coverage", type=float, default=DEFAULT_MIN_COVERAGE)
    ap.add_argument("--major-rhigh", type=float, help="m above MHW; default the median annual maximum")
    ap.add_argument("--major-hours", type=float, help="hours above the berm; default the median annual maximum")
    ap.add_argument("--variant", default=env.DEFAULT_STORM_VARIANT)
    a = ap.parse_args(argv)

    if a.window_dates:
        windows = [Window.from_dates(*a.window_dates)]
    elif a.window:
        windows = [Window.from_years(*a.window)]
    else:
        windows = [Window.centred_on(k) for k in (a.centred_on or sorted(SURVEY_ANCHORS))]

    fs.apply_style()
    events = load_events(a.variant)
    major = {"rhigh": a.major_rhigh if a.major_rhigh is not None
             else annual_max_median(events, "rhigh_m", sf.BERM_MHW),
             "hours": a.major_hours if a.major_hours is not None
             else annual_max_median(events, "raw_hours", 0)}
    print(f"{len(events)} events 1984-2024; major = Rhigh >= {major['rhigh']:.2f} m MHW "
          f"or >= {major['hours']:.0f} h above the berm")
    for w in windows:
        run(w, events, major, a.variant, a.context_years, a.recovery_days, a.min_coverage)


if __name__ == "__main__":
    main()
