"""
storm_construction_figures.py
==============================================================================
How the storm series the model reads is built, drawn from the same records by
the same rules as the generator
(scripts/input_prep/3-env-forcings/3-storms/historical_storm_creation_v3_HAT.py).

    python scripts/figure_making/pipeline/3-storms/storm_construction_figures.py

Writes to output/figures/3-model-inputs/3-forcing/:

    storm_construction_steps.png   the chain on a worked stretch (Edouard and
                                   Fran, late August - early September 1996):
                                   Duck water level and WIS waves -> Stockdon
                                   run-up -> total water level against the berm
                                   -> hours above it grouped into events ->
                                   events split where the storm hours break ->
                                   each event trimmed to 24 h around its peak
    storm_events_by_duration.png   every event of 1996-2024 by its length above
                                   the berm and its Rhigh, kept whole or
                                   trimmed; and the storm hours each year
                                   against the hours the model receives

THE RULE (the series in use, v3_split12_trim24, adopted 2026-09-29)
    1 hours when total water level exceeds the berm are storm hours
    2 storm hours less than 24 h apart are grouped into one event
    3 an event is split wherever consecutive storm hours are >= 12 h apart; a
      piece shorter than 8 h is folded into the piece before it (after, for
      the first)
    4 an event of fewer than 8 storm hours is not a storm
    5 an event longer than 24 h is cut to the 24 storm hours centred on its
      peak total water level; Rhigh, Rlow and the period come from what is kept

THE REPRODUCTION
    The generator is a script with module-level execution, so its logic is
    re-implemented here (build_events) and CHECKED before any figure is drawn:
    the events must equal, row for row, the committed
    <window>_storms_v3_split12_trim24_summary.csv of 1996_2010 and 2010_2024.
    The run stops if they do not.

REWORKED 2026-09-29 (Hannah: "rework the construction figures"). Until then
this drew the v3_72 rule -- events over 72 h dropped whole -- which has not
been the model's input since 2026-09-28 (trim24) and did not have the split.
"""

from __future__ import annotations

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
    apply_style, C, C_1997, INK, INK_MUTED, figsize, figure_dir, save, record_caption,
    _title, open_frame,
)

OUT = figure_dir("inputs", "3-forcing")

# the generator's inputs (historical_storm_creation_v3_HAT.py, "user inputs",
# and the command line of the series in use)
BEACH_SLOPE = 0.06
BERM_NAVD = 1.7
MHW_NAVD = 0.36
WEATHER_GROUPING_H = 24
SPLIT_GAP_H = 12
MIN_DUR_H = 8
TRIM_H = 24
VARIANT = "v3_split12_trim24"
WINDOWS = ((1996, 2010), (2010, 2024))
# the worked stretch: Edouard (peak 1 Sep) and Fran (peak 6 Sep) 1996
EXAMPLE = (pd.Timestamp("1996-08-27"), pd.Timestamp("1996-09-08"))

C_KEPT = C["ACCENT"]        # storm hours the model receives
C_CUT = C["BASE"]           # storm hours the trim removes
C_SHORT = C["BASE_FILL"]    # runs that never make an 8 h storm


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


def _split(group):
    """The generator's split: cut where consecutive storm hours are >=
    SPLIT_GAP_H apart, fold a piece under MIN_DUR_H into its neighbour."""
    group = group.sort_values("Time")
    gaps = group["Time"].diff().dt.total_seconds().div(3600).fillna(0).values
    cuts = [i for i, g in enumerate(gaps) if g >= SPLIT_GAP_H]
    bounds = [0] + cuts + [len(group)]
    pieces = [list(range(a, b)) for a, b in zip(bounds[:-1], bounds[1:])]
    i = 0
    while len(pieces) > 1 and i < len(pieces):
        if len(pieces[i]) < MIN_DUR_H:
            j = i - 1 if i > 0 else i + 1
            pieces[j] = sorted(pieces[j] + pieces[i])
            del pieces[i]
            i = 0
            continue
        i += 1
    return [group.iloc[p] for p in pieces]


def build_events(df, window_start_year):
    """create_storms() with the series-in-use options. Returns the hourly
    table (with the grouped system id) and EVERY piece, kept or not, with the
    hours the trim keeps. Units as the generator: Rhigh/Rlow in dam above MHW,
    duration in hours above the berm."""
    df = df.copy()
    df["R2"] = r2_stockdon(df.Hs, df.Tp)
    df["TWL"] = df.water_level + df.R2
    d = pd.DataFrame({"Time": pd.to_datetime(df.index), "TWL": df.TWL.values, "Tp": df.Tp.values})
    d = d.sort_values("Time").reset_index(drop=True)
    d["AboveBerm"] = d.TWL > BERM_NAVD
    d["StormStart"] = d.AboveBerm & (~d.AboveBerm.shift(1, fill_value=False))
    above = d.index[d.AboveBerm]
    for prev, cur in zip(above[:-1], above[1:]):
        if d.at[cur, "StormStart"]:
            if (d.Time[cur] - d.Time[prev]).total_seconds() / 3600 < WEATHER_GROUPING_H:
                d.at[cur, "StormStart"] = False
    d["StormID"] = d.StormStart.cumsum().astype(float)
    d.loc[~d.AboveBerm, "StormID"] = np.nan
    rows = []
    for sid, g in d.dropna(subset=["StormID"]).groupby("StormID"):
        for k, piece in enumerate(_split(g)):
            raw = len(piece)
            start = piece.Time.iloc[0]
            kept = piece
            if raw > TRIM_H:
                i = int(np.argmax(piece.TWL.values))
                lo = min(max(0, i - TRIM_H // 2), raw - TRIM_H)
                kept = piece.iloc[lo:lo + TRIM_H]
            peak = kept.TWL.idxmax()
            rows.append(dict(
                system=sid, piece=k, calendar_year=start.year,
                StartTime=kept.Time.iloc[0], EndTime=kept.Time.iloc[-1],
                raw_start=start, raw_end=piece.Time.iloc[-1],
                Rhigh=(kept.TWL.max() - MHW_NAVD) / 10, Rlow=(kept.TWL.min() - MHW_NAVD) / 10,
                period=kept.loc[peak, "Tp"], duration=len(kept), raw_hours=raw,
                trimmed_from=raw if raw > TRIM_H else 0, storm=raw >= MIN_DUR_H,
                kept_times=kept.Time.values))
    ev = pd.DataFrame(rows)
    ev["time"] = ev.calendar_year - int(window_start_year) + 1
    return d, ev


def verify():
    """The reproduction must equal the committed series row for row."""
    lines = []
    for w in WINDOWS:
        _, ev = build_events(load_merged(*w), w[0])
        k = ev[ev.storm].reset_index(drop=True)
        ref = pd.read_csv(env.storm_summary_file(*w, variant=VARIANT), parse_dates=["StartTime", "EndTime"])
        ok = (len(k) == len(ref)
              and (k.StartTime.values == ref.StartTime.values).all()
              and (k.EndTime.values == ref.EndTime.values).all()
              and np.allclose(k.Rhigh, ref.Rhigh) and np.allclose(k.Rlow, ref.Rlow)
              and np.allclose(k.period, ref.period) and (k.duration.values == ref.duration.values).all()
              and (k.trimmed_from.values == ref.trimmed_from.values).all()
              and (k.time.values == ref.time.values).all())
        lines.append(f"{w[0]}_{w[1]}: {len(k)} storms reproduced vs {len(ref)} in {VARIANT} -> "
                     f"{'IDENTICAL' if ok else 'DIFFERENT'}")
        if not ok:
            raise SystemExit("reproduction differs from the committed storm file:\n  " + "\n  ".join(lines))
    return lines


# =============================================================================
# FIGURES
# =============================================================================

def fig_steps():
    lo, hi = EXAMPLE
    df = load_merged(lo.year, lo.year)
    d, ev = build_events(df, lo.year)
    s = df.loc[lo:hi].copy()
    s["R2"] = r2_stockdon(s.Hs, s.Tp)
    s["TWL"] = s.water_level + s.R2
    dd = d[(d.Time >= lo) & (d.Time <= hi)]
    evw = ev[(ev.raw_end >= lo) & (ev.raw_start <= hi)]
    systems = evw.system.unique()

    fig, axes = plt.subplots(4, 1, figsize=figsize("double", height=7.4), sharex=True,
                             constrained_layout=True, gridspec_kw=dict(height_ratios=[1, 1, 1.3, 1.25]))
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
    ax.set_ylim(None, max(s.TWL.max() * 1.3, 3))
    open_frame(ax)
    _title(ax, 2, "Total water level = water level + run-up")

    # (d) three rows: grouped, split, trimmed
    ax = axes[3]
    y_group, y_split, y_trim = 2.0, 1.0, 0.0
    hours = dd[dd.AboveBerm]
    ax.plot(hours.Time, np.full(len(hours), y_group), "|", color=INK_MUTED, ms=9, mew=1.0)
    for sid in systems:
        g = hours[hours.StormID == sid]
        if len(g):
            ax.plot([g.Time.min(), g.Time.max()], [y_group + 0.28] * 2, color=INK_MUTED, lw=1.6,
                    solid_capstyle="butt")
    for _, e in evw.iterrows():
        col = C_KEPT if e.storm else C_SHORT
        ax.plot([e.raw_start, e.raw_end], [y_split] * 2, color=col, lw=4, solid_capstyle="butt")
        if e.storm:
            kt = pd.to_datetime(e.kept_times)
            piece_hours = hours[(hours.Time >= e.raw_start) & (hours.Time <= e.raw_end)]
            cut = piece_hours[~piece_hours.Time.isin(kt)]
            ax.plot(cut.Time, np.full(len(cut), y_trim), "|", color=C_CUT, ms=9, mew=1.0)
            ax.plot(kt, np.full(len(kt), y_trim), "|", color=C_KEPT, ms=9, mew=1.0)
            mid = kt.min() + (kt.max() - kt.min()) / 2
            lab = (f"{e.raw_hours} h -> {e.duration} h, Rhigh {e.Rhigh * 10:.2f} m"
                   if e.trimmed_from else f"{e.duration} h, Rhigh {e.Rhigh * 10:.2f} m")
            ax.text(mid, y_trim - 0.38, lab, ha="center", va="top", fontsize=7, color=INK)
    ax.set_yticks([y_trim, y_split, y_group])
    ax.set_yticklabels([f"3. trimmed to {TRIM_H} h", f"2. split at\n>= {SPLIT_GAP_H} h gaps",
                        f"1. grouped\n(< {WEATHER_GROUPING_H} h apart)"], fontsize=7)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    open_frame(ax)
    ax.set_ylim(-1.0, 3.25)
    ax.legend(handles=[Line2D([], [], color=C_KEPT, lw=3, label="storm hours the model receives"),
                       Line2D([], [], color=C_CUT, lw=3, label="storm hours the trim removes"),
                       Line2D([], [], color=C_SHORT, lw=3, label=f"under {MIN_DUR_H} h: folded or not a storm")],
              frameon=False, fontsize=7, loc="upper center", ncol=3)
    _title(ax, 3, "Storm hours made into model events")
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))
    ax.set_xlabel(str(lo.year))
    ax.set_xlim(lo, hi)

    storms = evw[evw.storm].sort_values("raw_start")
    names = ["Edouard", "Fran"] if len(storms) == 2 else [f"event {i + 1}" for i in range(len(storms))]
    parts = [f"{n} ({e.raw_hours} h above the berm, kept as the {e.duration} h around its "
             f"{e.Rhigh * 10:.2f} m MHW peak)" for n, (_, e) in zip(names, storms.iterrows())]
    out = save(fig, OUT / "storm_construction_steps.png")
    plt.close(fig)
    record_caption(out[0],
        "How the model's storm series is built, on Hurricanes Edouard and Fran, late August to early September "
        "1996. (a) Hourly water level at Duck (NOAA 8651370, m NAVD88). (b) Significant wave height at WIS "
        "station ST63228; the peak period Tp is used too. (c) The 2% exceedance run-up of Stockdon et al. (2006) "
        "on a 0.06 beach slope (grey) added to the water level gives the total water level (black); hours when "
        f"it stands above the {BERM_NAVD} m NAVD88 berm (purple fill) are storm hours. (d) The storm hours "
        "(ticks) become model events in three steps. 1: hours less than "
        f"{WEATHER_GROUPING_H} h apart are grouped (grey bar), and because the berm is topped at most high "
        "tides here, the two storms chain into one group. 2: the group is split wherever consecutive storm "
        f"hours are {SPLIT_GAP_H} h or more apart; a piece under {MIN_DUR_H} h joins its neighbour, and a group "
        f"under {MIN_DUR_H} h is not a storm (pale). 3: an event longer than {TRIM_H} h is cut to the {TRIM_H} "
        "storm hours centred on its peak (purple kept, grey removed). Here: " + "; ".join(parts) + ". Each kept "
        "event is one row of the storm file: Rhigh and Rlow are its highest and lowest total water level less "
        "MHW (0.36 m NAVD88), in decametres; period is Tp at the peak; duration is its kept hours. Without "
        "step 2 the trim kept only Edouard, and Fran was not in the model. The reproduction used here returns "
        f"the committed 1996-2010 and 2010-2024 {VARIANT} files row for row.")
    return out


def fig_events():
    evs = []
    for w in WINDOWS:
        _, ev = build_events(load_merged(*w), w[0])
        last = w[1] - 1 if w == WINDOWS[0] else w[1]        # the calendar years each run spends
        evs.append(ev[ev.storm & (ev.calendar_year <= last)])
    ev = pd.concat(evs).reset_index(drop=True)
    ev["Rhigh_m"] = ev.Rhigh * 10
    whole, cut = ev[ev.trimmed_from == 0], ev[ev.trimmed_from > 0]

    fig, axes = plt.subplots(1, 2, figsize=figsize("double", height=3.4), constrained_layout=True,
                             gridspec_kw=dict(width_ratios=[1.2, 1]))
    ax = axes[0]
    ax.scatter(whole.raw_hours, whole.Rhigh_m, s=10, color=C["BASE"], lw=0,
               label=f"kept whole ({len(whole)})")
    ax.scatter(cut.raw_hours, cut.Rhigh_m, s=14, color=C_KEPT, lw=0,
               label=f"trimmed to {TRIM_H} h ({len(cut)})")
    ax.axvline(TRIM_H, color=INK_MUTED, lw=0.8, ls=(0, (3, 2)))
    for k, (_, e) in enumerate(ev.nlargest(5, "Rhigh_m").iterrows()):
        ax.annotate(f"{e.raw_start:%b %Y}", xy=(e.raw_hours, e.Rhigh_m), xytext=(5, 6 if k % 2 == 0 else -7),
                    textcoords="offset points", fontsize=7, va="center")
    ax.set_xscale("log")
    ax.set_xticks([8, 24, 48, 72, 120])
    ax.set_xticklabels(["8", "24", "48", "72", "120"])
    ax.set_xlabel("storm hours above the berm, before trimming (log)")
    ax.set_ylabel("Rhigh (m MHW)")
    open_frame(ax)
    ax.legend(frameon=False, fontsize=7.5, loc="lower right")
    _title(ax, 0, "Every event 1996-2024")

    ax = axes[1]
    years = np.arange(1996, 2025)
    raw = ev.groupby("calendar_year").raw_hours.sum().reindex(years, fill_value=0)
    kept = ev.groupby("calendar_year").duration.sum().reindex(years, fill_value=0)
    ax.bar(years, raw, color=C["BASE_FILL"], width=0.8, label="storm hours")
    ax.bar(years, kept, color=C_KEPT, width=0.8, label="hours the model receives")
    ax.set_xlabel("year")
    ax.set_ylabel("hours above the berm")
    open_frame(ax)
    ax.legend(frameon=False, fontsize=7.5, loc="upper left", ncol=2)
    ax.set_ylim(0, raw.max() * 1.18)
    _title(ax, 1, "Storm hours each year")

    out = save(fig, OUT / "storm_events_by_duration.png")
    plt.close(fig)
    top = ev.nlargest(5, "Rhigh_m")
    record_caption(out[0],
        f"Every storm event in the model's series ({VARIANT}), 1996-2024: {len(ev)} events, 1996-2009 from the "
        "1996-2010 run and 2010-2024 from the 2010-2024 run. (a) Each event by its storm hours above the berm "
        f"before trimming and its Rhigh. The dashed line is the {TRIM_H} h trim: {len(cut)} longer events reach "
        f"the model as the {TRIM_H} hours around their peak, and none is dropped. The five highest are labelled ("
        + ", ".join(f"{e.raw_start:%b %Y} {e.Rhigh_m:.2f} m" for _, e in top.iterrows())
        + f"). (b) Storm hours each year (pale) against the hours the model receives after the trim (purple): "
        f"{int(kept.sum())} of {int(raw.sum())} hours, {100 * kept.sum() / raw.sum():.0f}%. The trim keeps each "
        "storm's peak and about two tidal cycles around it; Barrier3D applies an event's peak Rhigh for every "
        "hour it lasts, so untrimmed events would multiply overwash volume "
        "(experiments/storms-and-overwash/2026-09-28-trim-length-adopted). Events are split before trimming, "
        f"wherever consecutive storm hours are {SPLIT_GAP_H} h or more apart, so a second storm in a chain "
        "keeps its own peak.")
    return out


def main():
    apply_style()
    for line in verify():
        print("  reproduction", line)
    for f in (fig_steps, fig_events):
        print(f.__name__, "->", f()[0].relative_to(REPO))


if __name__ == "__main__":
    main()
