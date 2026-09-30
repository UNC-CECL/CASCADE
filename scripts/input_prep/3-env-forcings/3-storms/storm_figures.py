"""
Storm record figures for the canonical 1996 -> 2010 -> 2024 chain, coloured by HURDAT2 type.

    python scripts/input_prep/3-env-forcings/3-storms/storm_figures.py
    python scripts/input_prep/3-env-forcings/3-storms/storm_figures.py --variant v3_72

Writes the record and characteristics figures, one folder per window. Details: scripts/input_prep/3-env-forcings/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""
import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

_REPO = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(_REPO / "scripts"))

from site_layer import hat_env_forcings as env  # noqa: E402
from site_layer.hat_figure_style import (C, C_1984, C_1997, SMOOTH_RAMP,  # noqa: E402
                                         INK, INK_MUTED, apply_style, caption, figsize,
                                         open_frame, save, support_dir, _title)

# --- CONFIG ------------------------------------------------------------------
MHW_NAVD88 = 0.36          # generator's MHW: 0 m MHW = 0.36 m NAVD88
BERM_NAVD88 = 1.7          # generator's berm_elevation
BERM_MHW = BERM_NAVD88 - MHW_NAVD88
PERIODS = ((1996, 2010), (2010, 2024))
PERIOD_COLOURS = (C_1984, C_1997)          # earlier red, later blue (house pair)
OUT_DIR = env.HINDCAST_STORMS.parent / "figures" / "1996_2024"   # figures/ holds one folder per window

# STORM TYPE, FROM THE NHC BEST TRACKS (HURDAT2)
HURDAT2_FILE = (_REPO / "data" / "hatteras_init" / "3-env-forcings" / "1-records" / "hurdat2"
                / "hurdat2-1851-2025-091226.txt")
CAPE_HATTERAS = (35.25, -75.53)
TC_STATUSES = ("TD", "TS", "HU", "SD", "SS")
TC_NEAR_KM = 500
TC_WINDOW_H = 24
TYPES = ("tropical", "other")
TYPE_LABELS = {"tropical": f"Tropical cyclone within {TC_NEAR_KM} km",
               "other": "Other high-water event"}
# House colours (Hannah, 2026-09-29)
TYPE_COLOURS = {"tropical": C["ADDED"], "other": SMOOTH_RAMP[1]}

# The generator's total water level, rebuilt to find each event's peak hour
BEACH_SLOPE = 0.06
# -----------------------------------------------------------------------------


# Duck water level and WIS waves on one hourly index
def load_forcing():
    wl = pd.read_csv(env.DUCK_GAUGE_FILE, index_col="t", parse_dates=True)["v"]
    wis = pd.read_csv(env.WIS_FILE, index_col="time", parse_dates=True)
    m = pd.DataFrame({"wl": wl, "Hs": wis.waveHs, "Tp": wis.waveTp,
                      "dir": wis.waveMeanDirection}).dropna()
    sq = np.sqrt(m.Hs * 9.81 * m.Tp ** 2 / (2 * np.pi))
    m["TWL"] = m.wl + 1.1 * (0.35 * BEACH_SLOPE * sq
                             + np.sqrt((0.75 * BEACH_SLOPE * sq) ** 2 + (0.06 * sq) ** 2) / 2)
    return m


# Both windows' storm summaries as one record
def load_record(variant, forcing):
    parts = []
    for (a, b), keep_to in zip(PERIODS, (PERIODS[0][1] - 1, PERIODS[1][1])):
        df = pd.read_csv(env.storm_summary_file(a, b, variant), parse_dates=["StartTime", "EndTime"])
        df = df[df.calendar_year <= keep_to].copy()
        df["pi"] = PERIODS.index((a, b))   # which run spends the event
        parts.append(df)
    # "period" in the summary is the wave period (Tp, s), not the hindcast period
    df = pd.concat(parts, ignore_index=True).rename(columns={"period": "tp_s"})
    df["peak"] = [forcing.TWL.loc[s:e].idxmax() for s, e in zip(df.StartTime, df.EndTime)]
    df["rhigh_m"] = df.Rhigh * 10.0
    rebuilt = forcing.TWL.loc[df.peak].values - MHW_NAVD88
    if not np.allclose(rebuilt, df.rhigh_m, atol=0.01):
        raise RuntimeError("rebuilt TWL does not reproduce the summary's Rhigh; the peak hours are wrong")
    trimmed = df["trimmed_from"] if "trimmed_from" in df else pd.Series(0, index=df.index)
    df["raw_hours"] = np.where(trimmed > 0, trimmed, df.duration)
    df["frac_year"] = df.peak.dt.year + (df.peak.dt.dayofyear - 1 + df.peak.dt.hour / 24) / 365.25
    return df


# Tropical/subtropical fixes since `first_year`, with the distance (km) of each fix from Cape Hatteras
def load_hurdat2(path=HURDAT2_FILE, first_year=1990):
    rows, name = [], None
    with open(path) as f:
        for line in f:
            p = [x.strip() for x in line.split(",")]
            if p[0][:2].isalpha():                     # header: AL092011, IRENE, 39,
                name = p[1].title()
                continue
            if int(p[0][:4]) < first_year or p[3] not in TC_STATUSES:
                continue
            lat = float(p[4][:-1]) * (1 if p[4][-1] == "N" else -1)
            lon = float(p[5][:-1]) * (-1 if p[5][-1] == "W" else 1)
            rows.append((name, pd.Timestamp(p[0] + p[1]), p[3], lat, lon))
    h = pd.DataFrame(rows, columns=["tc_name", "t", "status", "lat", "lon"])
    la, lo = np.radians(h.lat), np.radians(h.lon)
    la0, lo0 = np.radians(CAPE_HATTERAS[0]), np.radians(CAPE_HATTERAS[1])
    h["km"] = 6371.0 * 2 * np.arcsin(np.sqrt(np.sin((la - la0) / 2) ** 2
                                             + np.cos(la) * np.cos(la0) * np.sin((lo - lo0) / 2) ** 2))
    return h


# Each event's type from HURDAT2: tropical near, tropical far, or nor'easter
def classify(df, forcing):
    h = load_hurdat2()
    win = pd.Timedelta(hours=TC_WINDOW_H)
    names, statuses, kms = [], [], []
    for t in df.peak:
        c = h[(h.t >= t - win) & (h.t <= t + win)]
        r = c.loc[c.km.idxmin()] if not c.empty else None
        names.append(r.tc_name if r is not None else None)
        statuses.append(r.status if r is not None else None)
        kms.append(r.km if r is not None else np.nan)
    df["tc_name"], df["tc_status"], df["tc_km"] = names, statuses, kms
    df["type"] = np.where(df.tc_km <= TC_NEAR_KM, "tropical", "other")
    df["wave_dir"] = forcing.dir.loc[df.peak].values
    df["name"] = np.where(df.type == "tropical", df.tc_name, None)
    return df


# A named storm's label
def label_text(row):
    if row["name"] == "Unnamed":    # HURDAT2's unnamed systems (the October 2000 subtropical storm)
        kind = "Subtropical" if str(row.get("tc_status", "")).startswith("S") else "Tropical"
        return f"{kind} storm {row.peak.year}"
    return f"{row['name']} {row.peak.year}"


# The line between the two windows, labelled
def period_boundary(ax, label_y=None):
    b = PERIODS[1][0]
    ax.axvline(b, color=INK_MUTED, lw=0.7, ls=(0, (4, 2)), zorder=1)
    if label_y is not None:
        for (a, e), x, ha in zip(PERIODS, (b - 0.3, b + 0.3), ("right", "left")):
            ax.text(x, label_y, f"{a}–{e} run", transform=ax.get_xaxis_transform(), ha=ha,
                    va="top", color=INK_MUTED, fontsize=7.5)


# Figure 1: the record

# Every event through time, and the count per year
def fig_record(df, variant, out_dir):
    y0, y1 = PERIODS[0][0], PERIODS[1][1]
    years = np.arange(y0, y1 + 1)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize("double", height=5.4), sharex=True,
                                   gridspec_kw=dict(height_ratios=[1, 2.6], hspace=0.1))
    halo = [pe.withStroke(linewidth=2.4, foreground="white")]

    # (a) events per year, other high-water events below, tropical above
    counts = {k: df[df.type == k].groupby(df.peak.dt.year).size().reindex(years, fill_value=0)
              for k in TYPES}
    bottom = np.zeros(len(years))
    for k in ("other", "tropical"):
        ax1.bar(years + 0.5, counts[k], width=0.74, bottom=bottom, color=TYPE_COLOURS[k],
                lw=0, zorder=2)
        bottom = bottom + counts[k].values
    total = sum(counts.values())
    means = []
    for i, (a, b) in enumerate(PERIODS):
        yrs = years[(years >= a) & (years < b + (i == 1))]
        m = total.loc[yrs].mean()
        means.append(m)
        ax1.plot([yrs[0] + 0.1, yrs[-1] + 0.9], [m, m], color=INK, lw=0.9, ls=(0, (5, 2)), zorder=3)
        # in the headroom above the bars, at the period's left edge
        ax1.text(yrs[0] + 0.3, 0.98, f"mean {m:.1f} events yr$^{{-1}}$ (dashed)",
                 transform=ax1.get_xaxis_transform(), ha="left", va="top", fontsize=7,
                 color=INK, path_effects=halo, zorder=4)
    ax1.set_ylabel("Events per year")
    ax1.set_ylim(0, total.max() * 1.3)
    ax1.yaxis.grid(True, zorder=0)
    ax1.set_axisbelow(True)
    open_frame(ax1)
    period_boundary(ax1)
    _title(ax1, 0, "")

    # (b) every event, coloured by type; tropical drawn last so they sit on top
    for k, size in (("other", 15), ("tropical", 19)):
        d = df[df.type == k]
        ax2.scatter(d.frac_year, d.rhigh_m, s=size, color=TYPE_COLOURS[k], alpha=0.9,
                    edgecolors="white", linewidths=0.35, zorder=3)
    ax2.axhline(BERM_MHW, color=INK_MUTED, lw=0.8, ls=(0, (4, 2)), zorder=1)
    ax2.text(y1 + 0.9, BERM_MHW + 0.04, f"berm crest ({BERM_MHW:.2f} m)", ha="right",
             fontsize=7, color=INK_MUTED, va="bottom", path_effects=halo)

    # EVERY tropical event is labelled (Hannah, 2026-09-29), once per storm
    top = (df[df["name"].notna()].sort_values("rhigh_m", ascending=False)
           .drop_duplicates(subset=["name", "calendar_year"]).sort_values("frac_year"))
    labelled = []
    for _, r in top.iterrows():
        txt = label_text(r)
        labelled.append(dict(peak=r.peak, rhigh_m_mhw=round(r.rhigh_m, 3),
                             rhigh_m_navd88=round(r.rhigh_m + MHW_NAVD88, 3),
                             hours_above_berm=int(r.raw_hours), label=txt,
                             type=r["type"], nearest_tc=r["tc_name"],
                             nearest_tc_km=round(r["tc_km"]) if np.isfinite(r["tc_km"]) else None))
    # The highest non-tropical event, marked by type (it has no official name)
    ne = df[df.type == "other"].nlargest(1, "rhigh_m").iloc[0]
    ax2.annotate("nor'easter", (ne.frac_year, ne.rhigh_m), xytext=(0, 5), textcoords="offset points",
                 ha="center", va="bottom", fontsize=6.8, fontstyle="italic", color=INK, zorder=6,
                 path_effects=[pe.withStroke(linewidth=2.2, foreground="white")])

    ax2.set_xlim(y0, y1 + 1)
    ax2.set_ylim(1.0, df.rhigh_m.max() + 0.45)
    ax2.set_ylabel(r"Peak total water level, $R_{high}$ (m above MHW)")
    ax2.set_xlabel("Year")
    ax2.set_xticks(np.arange(y0, y1 + 2, 2))
    ax2.set_xticks(np.arange(y0, y1 + 2, 1), minor=True)
    ax2.yaxis.grid(True, zorder=0)
    ax2.set_axisbelow(True)
    open_frame(ax2)
    period_boundary(ax2, label_y=0.99)
    _title(ax2, 1, "")
    _place_labels(ax2, top, [d["label"] for d in labelled], df)

    # one legend for both panels: the storm types
    handles = [Patch(color=TYPE_COLOURS[k], label=TYPE_LABELS[k]) for k in TYPES]
    fig.legend(handles=handles, loc="lower right", bbox_to_anchor=(0.9, 0.885), ncol=3,
               fontsize=7.5, handlelength=1.2)

    n_trim = int((df.raw_hours > 24).sum())
    caption(fig, (
        f"Storm record for the 1996–2024 hindcast: the {variant} series the model runs on "
        f"({len(df)} events; 1996–2009 from the 1996–2010 file, 2010–2024 from the 2010–2024 file, "
        f"the calendar years each run spends). Events are hours when total water level (Duck gauge "
        f"8651370 + Stockdon 2006 R2% from WIS 63228 waves, beach slope 0.06) exceeds the berm "
        f"({BERM_NAVD88} m NAVD88 = {BERM_MHW:.2f} m MHW), grouped when <24 h apart, kept when ≥8 h; "
        + (f"each group split where the water stays below the berm ≥{variant.split('split')[1].split('_')[0]} h, "
           if "split" in variant else "")
        + f"{n_trim} events longer than 24 h reach the model as the 24 h around their peak. "
        f"Dark orange: a tropical or subtropical cyclone in the NHC best tracks (HURDAT2) was within "
        f"{TC_NEAR_KM} km of Cape Hatteras within {TC_WINDOW_H} h of the event's peak "
        f"({(df.type == 'tropical').sum()} events). Slate: every other event "
        f"({(df.type == 'other').sum()}), mostly nor'easters "
        f"({int(df[df.type == 'other'].peak.dt.month.isin([10, 11, 12, 1, 2, 3, 4]).sum())} of them "
        f"peak October–April), but identified only by the absence of a nearby cyclone. Water levels "
        f"are the model's estimate from the gauge and hindcast waves, not observations at Hatteras. "
        f"Each event is placed at the hour of its peak total water level. "
        f"(a) Events per calendar year; dashed lines are each period's mean "
        f"({means[0]:.1f} yr⁻¹ 1996–2009, {means[1]:.1f} yr⁻¹ 2010–2024). "
        f"(b) Rhigh of every event (m above MHW; add {MHW_NAVD88} m for NAVD88). The dashed vertical "
        f"line is the 2010 boundary between the two hindcast runs. Every tropical cyclone is labelled with "
        f"its HURDAT2 name, at its highest event where it made two; the other events have no official "
        f"name and are unlabelled, except the highest, the nor'easter of {ne.peak:%B %Y}."
        + ("" if "split" in variant else " Fran 1996 and Jose 2017 are absent: each fell inside a longer "
           "grouped event (with Edouard and with Maria) whose 24 h trim kept the other storm's peak.")))
    save(fig, out_dir / "storm_record_1996_2024.png", bbox_inches="tight")
    plt.close(fig)
    return pd.DataFrame(labelled)


# Label placement: first free offset, highest events first
_LABEL_CANDIDATES = [(0, 5), (0, -5), (14, 6), (-14, 6), (14, -6), (-14, -6), (0, 14), (0, -14),
                     (24, 12), (-24, 12), (24, -12), (-24, -12), (0, 24), (0, -24), (34, 0), (-34, 0),
                     (30, 22), (-30, 22), (30, -22), (-30, -22), (0, 34), (0, -34)]
LABEL_PT = 6.5


# Place event labels clear of each other and the dots
def _place_labels(ax, top, texts, events):
    from matplotlib.text import Text
    fig = ax.figure
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    ax_box = ax.get_window_extent(rend)
    dots = ax.transData.transform(np.c_[events.frac_year, events.rhigh_m])
    r_dot = 3.2 * fig.dpi / 72                         # a dot's radius, display units
    placed = [t.get_window_extent(rend) for t in ax.texts]  # period, berm and nor'easter text
    # the berm line and the 2010 divider, as thin bands a label must not sit on
    from matplotlib.transforms import Bbox
    (bx, by), = ax.transData.transform([(PERIODS[1][0], BERM_MHW)])
    pad = 2.0 * fig.dpi / 72
    placed += [Bbox([[ax_box.x0, by - pad], [ax_box.x1, by + pad]]),
               Bbox([[bx - pad, ax_box.y0], [bx + pad, ax_box.y1]])]
    order = np.argsort(-top.rhigh_m.values)
    rows = list(top.iterrows())
    for i in order:
        (_, r), txt = rows[i], texts[i]
        own = ax.transData.transform((r.frac_year, r.rhigh_m))
        best = None
        for dx, dy in _LABEL_CANDIDATES:
            ann = ax.annotate(txt, (r.frac_year, r.rhigh_m), xytext=(dx, dy), textcoords="offset points",
                              ha="center", va="bottom" if dy >= 0 else "top", fontsize=LABEL_PT)
            ann.update_positions(rend)     # an annotation is only positioned at draw time
            bb = Text.get_window_extent(ann, rend).expanded(1.06, 1.2)
            ann.remove()
            inside = ax_box.x0 <= bb.x0 and bb.x1 <= ax_box.x1 and ax_box.y0 <= bb.y0 and bb.y1 <= ax_box.y1
            hit_lab = sum(bb.overlaps(o) for o in placed)
            near = ((dots[:, 0] > bb.x0 - r_dot) & (dots[:, 0] < bb.x1 + r_dot)
                    & (dots[:, 1] > bb.y0 - r_dot) & (dots[:, 1] < bb.y1 + r_dot))
            near &= np.hypot(*(dots - own).T) > 1.0          # not its own dot
            cost = (0 if inside else 1e6) + 1e4 * hit_lab + 6 * int(near.sum()) + np.hypot(dx, dy)
            if best is None or cost < best[0]:
                best = (cost, dx, dy, bb)
            if cost < 6:
                break
        _, dx, dy, bb = best
        placed.append(bb)
        ax.annotate(txt, (r.frac_year, r.rhigh_m), xytext=(dx, dy), textcoords="offset points",
                    ha="center", va="bottom" if dy >= 0 else "top", fontsize=LABEL_PT, color=INK, zorder=6,
                    arrowprops=(dict(arrowstyle="-", lw=0.4, color=INK_MUTED, shrinkA=0, shrinkB=2.5)
                                if np.hypot(dx, dy) > 10 else None),
                    path_effects=[pe.withStroke(linewidth=2.0, foreground="white")])


# Figure 2: the two periods compared

# Height, duration and seasonality by type
def fig_characteristics(df, variant, out_dir):
    fig, axs = plt.subplots(2, 2, figsize=figsize("double", height=5.4),
                            gridspec_kw=dict(hspace=0.45, wspace=0.3))
    (a, b), (c, d) = axs
    nyears = [PERIODS[0][1] - PERIODS[0][0], PERIODS[1][1] - PERIODS[1][0] + 1]
    labels = [f"{p[0]}–{p[1]} ({n} yr)" for p, n in zip(PERIODS, nyears)]
    labels = [f"{PERIODS[0][0]}–{PERIODS[0][1] - 1} ({nyears[0]} yr, 1996–2010 run)",
              f"{PERIODS[1][0]}–{PERIODS[1][1]} ({nyears[1]} yr, 2010–2024 run)"]

    # (a) exceedance: events per year with Rhigh >= x
    x = np.linspace(BERM_MHW, df.rhigh_m.max() + 0.1, 300)
    for i in (0, 1):
        r = df[df.pi == i].rhigh_m.values
        a.step(x, [(r >= v).sum() / nyears[i] for v in x], where="post", color=PERIOD_COLOURS[i],
               lw=1.3, label=labels[i])
    a.set_yscale("log")
    a.set_xlabel(r"$R_{high}$ (m MHW)")
    a.set_ylabel(r"Events per year with $R_{high}$ ≥ x")
    a.set_xlim(BERM_MHW, df.rhigh_m.max() + 0.1)
    a.yaxis.grid(True, which="major")
    a.xaxis.grid(True)
    open_frame(a)
    _title(a, 0, "")

    # (b) seasonality: events per year by month
    months = np.arange(1, 13)
    wbar = 0.38
    for i in (0, 1):
        cnt = df[df.pi == i].groupby(df.peak.dt.month).size().reindex(months, fill_value=0) / nyears[i]
        b.bar(months + (i - 0.5) * wbar, cnt, width=wbar, color=PERIOD_COLOURS[i], zorder=2)
    b.set_xticks(months)
    b.set_xticklabels(list("JFMAMJJASOND"))
    b.set_xlabel("Month of peak")
    b.set_ylabel("Events per year")
    b.yaxis.grid(True, zorder=0)
    b.set_axisbelow(True)
    open_frame(b)
    _title(b, 1, "")

    # (c) hours above the berm, before trimming
    bins = np.concatenate([np.arange(8, 49, 4), [72, 96, 144, 200]])
    for i in (0, 1):
        h = np.clip(df[df.pi == i].raw_hours.values, None, bins[-1] - 1)
        cnt, _ = np.histogram(h, bins=bins)
        c.stairs(cnt / nyears[i], np.arange(len(bins)), color=PERIOD_COLOURS[i], lw=1.3)
    c.set_xticks(np.arange(len(bins)))
    c.set_xticklabels([str(v) if v not in (200,) else "" for v in bins], fontsize=7)
    c.set_xlim(0, len(bins) - 1)
    for hrs, txt in ((24, "longer events\ntrimmed to 24 h"),):
        xi = list(bins).index(hrs)
        c.axvline(xi, color=INK_MUTED, lw=0.7, ls=(0, (4, 2)))
        c.text(xi + 0.12, 0.97, txt, transform=c.get_xaxis_transform(), fontsize=6.8,
               color=INK_MUTED, va="top")
    c.set_xlabel("Hours above the berm (before trimming; bins unequal)")
    c.set_ylabel("Events per year")
    c.yaxis.grid(True)
    open_frame(c)
    _title(c, 2, "")

    # (d) Rhigh vs peak wave period
    for i in (0, 1):
        e = df[df.pi == i]
        d.scatter(e.tp_s, e.rhigh_m, s=10, color=PERIOD_COLOURS[i], alpha=0.7,
                  edgecolors="white", linewidths=0.3)
    d.axhline(BERM_MHW, color=INK_MUTED, lw=0.7, ls=(0, (4, 2)))
    d.set_xlabel("Peak wave period, $T_p$ (s)")
    d.set_ylabel(r"$R_{high}$ (m MHW)")
    d.yaxis.grid(True)
    d.xaxis.grid(True)
    open_frame(d)
    _title(d, 3, "")

    fig.legend(handles=[Line2D([], [], color=PERIOD_COLOURS[i], lw=1.6, label=labels[i]) for i in (0, 1)],
               loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.02), fontsize=8)
    fig.subplots_adjust(bottom=0.14)

    rates = [len(df[df.pi == i]) / nyears[i] for i in (0, 1)]
    caption(fig, (
        f"The two hindcast periods' storm series compared ({variant}; 1996–2009 as spent by the "
        f"1996–2010 run, red, {rates[0]:.1f} events/yr; 2010–2024 as spent by the 2010–2024 run, blue, "
        f"{rates[1]:.1f} events/yr). Every panel is normalised per year so the 14- and 15-year periods compare. "
        f"(a) Exceedance: events per year whose Rhigh reaches the value on the x-axis, from the berm "
        f"({BERM_MHW:.2f} m MHW) up; log scale. "
        f"(b) Events per year by the month of the peak. "
        f"(c) Length of each event above the berm as it came out of grouping, before the 24 h trim "
        f"(bins widen past 48 h; the last bin holds everything ≥144 h). Events right of the dashed line "
        f"reach the model as the 24 h around their peak. "
        f"(d) Rhigh against the event's wave period (Tp, the value the model receives)."))
    save(fig, out_dir / "storm_characteristics_1996_2024.png", bbox_inches="tight")
    plt.close(fig)


# Run: load, classify, draw both figures
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--variant", default=env.DEFAULT_STORM_VARIANT)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()

    apply_style()
    forcing = load_forcing()
    df = classify(load_record(args.variant, forcing), forcing)
    labelled = fig_record(df, args.variant, args.out_dir)
    fig_characteristics(df, args.variant, args.out_dir)
    labelled.to_csv(support_dir(args.out_dir) / "storm_record_1996_2024_labelled_events.csv", index=False)
    types = df[["peak", "calendar_year", "rhigh_m", "raw_hours", "wave_dir", "type",
                "tc_name", "tc_status", "tc_km"]].copy()
    types["tc_km"] = types.tc_km.round(0)
    types.round(3).to_csv(support_dir(args.out_dir) / "storm_types_1996_2024.csv", index=False)
    print(df.type.value_counts().to_string())
    print(f"{len(df)} events; figures in {args.out_dir}")
    print(labelled.to_string(index=False))


if __name__ == "__main__":
    main()
