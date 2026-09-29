"""
storm_construction_figures.py
==============================================================================
How the storm series the model reads is built, drawn from the same records by
the same rules as the generator
(scripts/input_prep/3-env-forcings/3-storms/historical_storm_creation_v3_HAT.py).

    python scripts/figure_making/pipeline/3-storms/storm_construction_figures.py [--max-dur 72]

Writes to output/figures/3-model-inputs/3-forcing/:

    storm_construction_steps.png   the chain on a worked month (September 2003):
                                   Duck water level and WIS waves -> Stockdon
                                   run-up -> total water level against the berm
                                   -> hours above it grouped into events ->
                                   events kept or dropped on duration
    storm_events_by_duration.png   every grouped event of 1996-2024 by duration
                                   and Rhigh, with the max-duration cut marked

THE REPRODUCTION
    The generator is a script with module-level execution, so its logic is
    re-implemented here line for line (build_events) and CHECKED before any
    figure is drawn: with --max-dur 72 the kept events must equal, row for
    row, the committed <window>_storms_v3_72_summary.csv of 1996_2010 and
    2010_2024. The run stops if they do not.

--max-dur sets the cut, so the figures can be redrawn for another variant.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

REPO = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))

from site_layer import hat_env_forcings as env  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    apply_style, C, C_1984, C_1997, INK, INK_MUTED, figsize, figure_dir, save, record_caption,
    _title, open_frame,
)

OUT = figure_dir("inputs", "3-forcing")

# the generator's inputs (historical_storm_creation_v3_HAT.py, "user inputs")
BEACH_SLOPE = 0.06
BERM_NAVD = 1.7
MHW_NAVD = 0.36
WEATHER_GROUPING_H = 24
MIN_DUR_H = 8
MAX_DUR_H = 72          # the committed series; --max-dur overrides
WINDOWS = ((1996, 2010), (2010, 2024))


# =============================================================================
# THE GENERATOR, REPRODUCED
# =============================================================================

def load_merged(start_year, end_year):
    """load_data(): an hourly index over the window, the Duck water level and
    WIS Hs/Tp placed on it, rows with any gap dropped."""
    idx = pd.date_range(f"{start_year}-01-01 00:00:00", f"{end_year}-12-31 23:00:00", freq="h")
    wl = pd.read_csv(env.DUCK_GAUGE_FILE, index_col="t")
    wl.index = pd.to_datetime(wl.index)
    wis = pd.read_csv(env.WIS_FILE, index_col="time")
    wis.index = pd.to_datetime(wis.index)
    df = pd.DataFrame(index=idx)
    df["water_level"] = wl["v"]
    df["Hs"] = wis["waveHs"]
    df["Tp"] = wis["waveTp"]
    return df.dropna()


def r2_stockdon(hs, tp, slope=BEACH_SLOPE):
    """calculate_r2_percent(): Stockdon et al. (2006) R2%."""
    L0 = 9.81 * tp ** 2 / (2 * np.pi)
    a = np.sqrt(hs * L0)
    setup = 0.35 * slope * a
    s_inc, s_ig = 0.75 * slope * a, 0.06 * a
    return 1.1 * (setup + np.sqrt(s_inc ** 2 + s_ig ** 2) / 2)


def build_events(df, max_dur=MAX_DUR_H, window_start_year=None):
    """create_storms(), every grouped event returned with `kept` and a reason,
    rather than only the kept ones. Units as the generator: Rhigh/Rlow in dam
    above MHW, duration in hours above the berm."""
    df = df.copy()
    df["R2"] = r2_stockdon(df.Hs, df.Tp)
    df["TWL"] = df.water_level + df.R2
    d = pd.DataFrame({"Time": pd.to_datetime(df.index), "TWL": df.TWL.values, "Tp": df.Tp.values})
    d = d.sort_values("Time").reset_index(drop=True)
    d["AboveBerm"] = d.TWL > BERM_NAVD
    d["StormStart"] = d.AboveBerm & (~d.AboveBerm.shift(1, fill_value=False))
    # continuous weather system: a run starting < 24 h after the previous
    # above-berm hour joins the previous storm
    above = d.index[d.AboveBerm]
    times = d.Time
    for prev, cur in zip(above[:-1], above[1:]):
        if d.at[cur, "StormStart"]:
            gap_h = (times[cur] - times[prev]).total_seconds() / 3600
            if gap_h < WEATHER_GROUPING_H:
                d.at[cur, "StormStart"] = False
    d["StormID"] = d.StormStart.cumsum().astype(float)
    d.loc[~d.AboveBerm, "StormID"] = np.nan
    rows = []
    for sid, g in d.dropna(subset=["StormID"]).groupby("StormID"):
        g = g.sort_values("Time")
        dur = len(g)
        peak = g.TWL.idxmax()
        rows.append(dict(
            calendar_year=g.Time.iloc[0].year, StartTime=g.Time.iloc[0], EndTime=g.Time.iloc[-1],
            Rhigh=(g.TWL.max() - MHW_NAVD) / 10, Rlow=(g.TWL.min() - MHW_NAVD) / 10,
            period=g.loc[peak, "Tp"], duration=dur,
            kept=MIN_DUR_H <= dur <= max_dur,
            reason=("kept" if MIN_DUR_H <= dur <= max_dur
                    else f"under {MIN_DUR_H} h" if dur < MIN_DUR_H else f"over {max_dur} h")))
    ev = pd.DataFrame(rows)
    kept = ev[ev.kept]
    anchor = int(kept.calendar_year.min()) if window_start_year is None else int(window_start_year)
    ev["time"] = ev.calendar_year - anchor + 1
    return d, ev


def verify(max_dur):
    """The reproduction must equal the committed 72 h files row for row."""
    lines = []
    for w in WINDOWS:
        df = load_merged(*w)
        _, ev = build_events(df, max_dur=72)
        k = ev[ev.kept].reset_index(drop=True)
        ref = pd.read_csv(env.storm_window_dir(*w) / f"{w[0]}_{w[1]}_storms_v3_72_summary.csv",
                          parse_dates=["StartTime", "EndTime"])
        ok = (len(k) == len(ref)
              and (k.StartTime.values == ref.StartTime.values).all()
              and (k.EndTime.values == ref.EndTime.values).all()
              and np.allclose(k.Rhigh, ref.Rhigh) and np.allclose(k.Rlow, ref.Rlow)
              and np.allclose(k.period, ref.period) and (k.duration.values == ref.duration.values).all()
              and (k.time.values == ref.time.values).all())
        lines.append(f"{w[0]}_{w[1]}: {len(k)} storms reproduced vs {len(ref)} in the file -> "
                     f"{'IDENTICAL' if ok else 'DIFFERENT'}")
        if not ok:
            raise SystemExit("reproduction differs from the committed storm file:\n  " + "\n  ".join(lines))
    return lines


# =============================================================================
# FIGURES
# =============================================================================

def fig_steps(max_dur):
    df = load_merged(2003, 2003)
    d, ev = build_events(df, max_dur=max_dur)
    lo, hi = pd.Timestamp("2003-08-30"), pd.Timestamp("2003-09-23")
    s = df.loc[lo:hi].copy()
    s["R2"] = r2_stockdon(s.Hs, s.Tp)
    s["TWL"] = s.water_level + s.R2
    dd = d[(d.Time >= lo) & (d.Time <= hi)]
    evw = ev[(ev.EndTime >= lo) & (ev.StartTime <= hi)]

    fig, axes = plt.subplots(4, 1, figsize=figsize("double", height=7.2), sharex=True,
                             constrained_layout=True, gridspec_kw=dict(height_ratios=[1, 1, 1.3, 0.9]))
    ax = axes[0]
    ax.plot(s.index, s.water_level, color=C_1997, lw=0.9)
    ax.set_ylabel("water level\n(m NAVD88)")
    open_frame(ax)
    _title(ax, 0, "Duck gauge (NOAA 8651370), hourly")
    ax = axes[1]
    ax.plot(s.index, s.Hs, color=INK, lw=0.9)
    ax.set_ylabel("H$_s$ (m)")
    open_frame(ax)
    _title(ax, 1, "WIS ST63228 significant wave height")
    ax = axes[2]
    ax.fill_between(s.index, s.water_level, s.TWL, color=C["BASE_FILL"], lw=0, label="run-up R2% (Stockdon 2006)")
    ax.plot(s.index, s.water_level, color=C_1997, lw=0.7, label="water level")
    ax.plot(s.index, s.TWL, color=INK, lw=0.9, label="total water level")
    ax.axhline(BERM_NAVD, color=C["ADDED"], lw=1.0, ls=(0, (3, 2)), label=f"berm {BERM_NAVD} m NAVD88")
    ax.fill_between(s.index, BERM_NAVD, s.TWL, where=s.TWL > BERM_NAVD, color=C["ACCENT_FILL"], lw=0,
                    label="hours above the berm")
    ax.set_ylabel("elevation\n(m NAVD88)")
    ax.legend(frameon=False, fontsize=7, ncol=3, loc="upper left")
    ax.set_ylim(None, max(s.TWL.max() * 1.25, 3))
    open_frame(ax)
    _title(ax, 2, "Total water level = water level + run-up")
    ax = axes[3]
    ids = dd.StormID.dropna().unique()
    for sid in ids:
        g = dd[dd.StormID == sid]
        row = ev.iloc[int(sid) - 1] if int(sid) - 1 < len(ev) else None
        col = (C["ACCENT"] if (row is not None and row.kept)
               else C["BASE"] if (row is not None and row.duration < MIN_DUR_H) else C_1984)
        ax.plot(g.Time, np.zeros(len(g)), "|", color=col, ms=14, mew=1.1)
    for _, e in evw[evw.duration >= MIN_DUR_H].iterrows():
        col = C["ACCENT"] if e.kept else C_1984
        ax.annotate("", xy=(e.EndTime, 0.55), xytext=(e.StartTime, 0.55),
                    arrowprops=dict(arrowstyle="-", color=col, lw=2.2))
        mid = e.StartTime + (e.EndTime - e.StartTime) / 2
        lab = (f"{e.duration} h above berm, Rhigh {e.Rhigh * 10:.2f} m MHW\n"
               + ("kept" if e.kept else f"DROPPED ({e.reason})"))
        ax.text(mid, 0.95, lab, ha="center", va="bottom", fontsize=7, color=INK)
    ax.set_yticks([])
    for sp in ("left",):
        ax.spines[sp].set_visible(False)
    open_frame(ax)
    ax.legend(handles=[Line2D([], [], color=C["ACCENT"], lw=2.2, label="event kept"),
                       Line2D([], [], color=C_1984, lw=2.2, label="event dropped (too long)"),
                       Line2D([], [], color=C["BASE"], lw=2.2, label=f"under {MIN_DUR_H} h: not a storm")],
              frameon=False, fontsize=7, loc="upper left", ncol=3)
    ax.set_ylim(-0.6, 2.8)
    _title(ax, 3, f"Hours above the berm, grouped when < {WEATHER_GROUPING_H} h apart")
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=3))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))
    ax.set_xlabel("2003")
    ax.set_xlim(lo, hi)

    isabel = evw.loc[evw.duration.idxmax()]
    out = save(fig, OUT / "storm_construction_steps.png")
    plt.close(fig)
    record_caption(out[0],
        "How the model's storm series is built, on September 2003. (a) Hourly water level at Duck (NOAA "
        "8651370, m NAVD88). (b) Significant wave height at WIS station ST63228; the peak period Tp is used "
        "too. (c) The 2% exceedance run-up of Stockdon et al. (2006) on a 0.06 beach slope (grey) added to the "
        f"water level gives the total water level (black); hours when it stands above the {BERM_NAVD} m NAVD88 "
        "berm (purple fill) are storm hours. (d) Runs of storm hours (bars) are chained into one event when they start "
        f"less than {WEATHER_GROUPING_H} h after the previous storm hour; an event is kept if it has "
        f"{MIN_DUR_H} to {max_dur} storm hours. Each kept event becomes one row of the storm file: Rhigh and "
        "Rlow are its highest and lowest total water level less MHW (0.36 m NAVD88), in decametres; period is "
        "Tp at the peak; duration is its hours above the berm. Hurricane Isabel's surge (18 September) chains "
        f"with the swell ahead of it into one {isabel.duration} h event reaching Rhigh "
        f"{isabel.Rhigh * 10:.2f} m MHW, which the {max_dur} h limit "
        f"{'drops entirely: the storm file has no Isabel' if not isabel.kept else 'keeps'}. The reproduction "
        "used here returns the committed 1996-2010 and 2010-2024 storm files row for row at a 72 h limit.")
    return out


def fig_events(max_dur):
    evs = []
    for w in WINDOWS:
        _, ev = build_events(load_merged(*w), max_dur=max_dur)
        evs.append(ev)
    ev = pd.concat(evs).drop_duplicates(subset=["StartTime"]).reset_index(drop=True)
    ev["Rhigh_m"] = ev.Rhigh * 10
    ev = ev[ev.duration >= MIN_DUR_H]
    kept, drop = ev[ev.kept], ev[~ev.kept]

    fig, axes = plt.subplots(1, 2, figsize=figsize("double", height=3.4), constrained_layout=True,
                             gridspec_kw=dict(width_ratios=[1.2, 1]))
    ax = axes[0]
    ax.scatter(kept.duration, kept.Rhigh_m, s=10, color=C["BASE"], lw=0, label=f"kept ({len(kept)})")
    ax.scatter(drop.duration, drop.Rhigh_m, s=18, color=C_1984, lw=0, label=f"dropped ({len(drop)})")
    ax.axvline(max_dur, color=INK_MUTED, lw=0.8, ls=(0, (3, 2)))
    for k, (_, e) in enumerate(drop.nlargest(5, "Rhigh_m").iterrows()):
        ax.annotate(f"{e.StartTime:%b %Y}", xy=(e.duration, e.Rhigh_m), xytext=(5, 6 if k % 2 == 0 else -7),
                    textcoords="offset points", fontsize=7, va="center")
    ax.set_xscale("log")
    ax.set_xticks([8, 24, 48, 72, 120, 200])
    ax.set_xticklabels(["8", "24", "48", "72", "120", "200"])
    ax.set_xlabel("storm hours above the berm (log)")
    ax.set_ylabel("Rhigh (m MHW)")
    open_frame(ax)
    ax.legend(frameon=False, fontsize=7.5, loc="lower right")
    _title(ax, 0, f"Every event 1996-2024; limit {max_dur} h")

    ax = axes[1]
    years = np.arange(1996, 2025)
    mk = kept.groupby("calendar_year").Rhigh_m.max().reindex(years)
    md = drop.groupby("calendar_year").Rhigh_m.max().reindex(years)
    ax.bar(years, mk, color=C["BASE"], width=0.8, label="largest kept")
    ax.scatter(years, md, color=C_1984, s=16, zorder=3, label="largest dropped")
    ax.set_xlabel("year")
    ax.set_ylabel("Rhigh (m MHW)")
    open_frame(ax)
    ax.set_ylim(0, 5.9)
    ax.legend(frameon=False, fontsize=7.5, loc="upper right", ncol=2)
    _title(ax, 1, "The largest storm each year")

    out = save(fig, OUT / "storm_events_by_duration.png")
    plt.close(fig)
    top = drop.nlargest(5, "Rhigh_m")
    record_caption(out[0],
        f"Every storm event the generator finds at Duck, 1996-2024, and which ones reach the model at a "
        f"{max_dur} h maximum duration. (a) Each event by its storm hours above the berm and its Rhigh; the "
        f"dashed line is the limit. {len(drop)} events are longer and are dropped whole, not shortened; the five "
        "highest are labelled (" + ", ".join(f"{e.StartTime:%b %Y} {e.Rhigh_m:.2f} m" for _, e in top.iterrows())
        + f"). The kept maximum is {kept.Rhigh_m.max():.2f} m. (b) The largest kept storm of each year (bars) "
        "against the largest dropped one (points): in the years with a point above its bar, the storm file "
        "is missing that year's largest event. Events from both hindcast windows' generator runs, 2010 "
        "counted once. Grouped runs under 8 h are not storms and are not shown.")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-dur", type=int, default=MAX_DUR_H)
    a = ap.parse_args()
    apply_style()
    for line in verify(a.max_dur):
        print("  reproduction", line)
    for f in (fig_steps, fig_events):
        print(f.__name__, "->", f(a.max_dur)[0].relative_to(REPO))


if __name__ == "__main__":
    main()
