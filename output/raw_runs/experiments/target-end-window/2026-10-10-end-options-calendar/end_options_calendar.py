#!/usr/bin/env python3
"""
The selected end window (last modelled year ±6 months) beside four other calendar-based ends.

    python end_options_calendar.py  ->  <period>/domain_net_change.csv, scores.csv, window_dates.csv,
                                        <period>/figures/end_options_<period>_<preset>.png,
                                        figures/end_options_both_periods_<preset>.png
                                        for <period> = calibration, test

No end is tied to a DEM date. Each option is defined from the last modelled year Y (calibration
2008, test 2025) and the model's final state (1 Jan Y+1); the same rule in both periods:
  selected    Y ±6 months (1 Jul Y-1 to 30 Jun Y+1), the advisor's pick
  year        Y only (1 Jan to 31 Dec Y)
  last_two    the last two modelled years (1 Jan Y-1 to 31 Dec Y)
  three_yr    Y ±1 yr, a 3-yr window (1 Jan Y-1 to 31 Dec Y+1)
  end_state   the model's final state ±1 yr (1 Jan Y to 1 Jan Y+2)
CoastSat ends 2026-01-13, so the test windows reaching past 2025 hold only two weeks of 2026.
Starts stay DEM-centred (the island offset's window). Interior GIS 2-89, both sides LOWESS-7.
No new runs.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-10
"""
from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = next(p for p in HERE.parents if (p / "pyproject.toml").exists())
NC = REPO / "scripts" / "input_prep" / "5-scr" / "3-rates" / "coastsat" / "net_change"
MS = REPO / "scripts" / "input_prep" / "5-scr" / "1-observations" / "mean_shoreline"
sys.path[:0] = [str(REPO / "scripts"), str(NC), str(MS)]
from coastsat_net_change import domain_net_change, smooth_like_model, window_means  # noqa: E402
from coastsat_mean_shoreline import timeseries_file  # noqa: E402
from coastsat_lrr import filter_dates, load_timeseries  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    DOMAIN_AXIS_LABEL, INK_MUTED, apply_style, figsize, mark_offaxis, offaxis_clause,
    open_frame, record_caption, save)
from cascade_pipeline.annotations import add_geographic_annotations  # noqa: E402
from site_layer.hat_observed_rates import NET_CHANGE_WINDOWS  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS, HATTERAS_DOMAINS as D  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
RAW = REPO / "output" / "raw_runs" / "matrix"
GIS = pd.RangeIndex(1, 91)
INTERIOR = slice(2, 89)
POS_HALF = 100.0                                # m, the net-change figures' fixed axis
# the same runs from the ±6-month island offset (offset v3), same tree under the 10-09 study
SIX_MO = HERE.parent / "2026-10-09-6mo-starts"
SIX_MO_INK, SIX_MO_DASH = "0.55", (0, (2, 1.2))
BOTH_HALF = 55.0                               # m, the two-period figure only, zoomed to show the options apart
DATA_END = pd.Timestamp("2026-01-13")           # last CoastSat position
SELECTED = "selected"


# The five ends, from the last modelled year Y: key -> (dates, legend name, colour, line width)
def ends_for(y):
    return {
        "selected": ((f"{y-1}-07-01", f"{y+1}-06-30"), f"Selected: {y} ±6 months", "#b2182b", 1.2),
        "year": ((f"{y}-01-01", f"{y}-12-31"), f"Last modelled year only ({y})", "#1b7837", 1.2),
        "last_two": ((f"{y-1}-01-01", f"{y}-12-31"), f"Last two modelled years ({y-1}-{y})", "#35978f", 1.2),
        "three_yr": ((f"{y-1}-01-01", f"{y+1}-12-31"), f"{y} ±1 yr, 3-yr window ({y-1}-{y+1})", "#762a83", 1.2),
        "end_state": ((f"{y}-01-01", f"{y+2}-01-01"),
                      f"Model end state (1 Jan {y+1}) ±1 yr", "#e08214", 1.2),
    }


PERIODS = {
    "calibration": dict(
        title="Calibration 1996-2009", start=NET_CHANGE_WINDOWS[(1996, 2009)][0],
        last_year=2008, model_start="1 Jan 1996",
        runs={
            "zeroBE, no groin": RAW / "1996_2009/zeroBE/HAT_1996_2009_zeroBE_offsetmetres_road_bdm_nogroin",
            "edgeBE, no groin": RAW / "1996_2009/edgeBE/HAT_1996_2009_edgeBE_offsetmetres_road_bdm_nogroin",
            "domainBE + groin": RAW / "1996_2009/domainBE/HAT_1996_2009_domainBE_offsetmetres_road_bdm_groinblock",
        }),
    "test": dict(
        title="Test 2009-2025", start=NET_CHANGE_WINDOWS[(2009, 2025)][0],
        last_year=2025, model_start="1 Jan 2009",
        runs={
            "zeroBE, no groin": RAW / "2009_2025/zeroBE/HAT_2009_2025_zeroBE_offsetmetres_road_bdm_nourish_nogroin",
            "edgeBE, no groin": RAW / "2009_2025/edgeBE/HAT_2009_2025_edgeBE_offsetmetres_road_bdm_nourish_nogroin",
            "domainBE + groin": RAW / "2009_2025/domainBE/HAT_2009_2025_domainBE_offsetmetres_road_bdm_nourish_groinblock",
        }),
}
# -----------------------------------------------------------------------------


def observed(start_w, end_w):
    s, e = window_means(start_w), window_means(end_w)
    both = s.index.intersection(e.index)
    tr = pd.DataFrame({"domain_number": s.loc[both, "domain_number"].astype(int),
                       "net_change_m": e.loc[both, "mean_chainage_m"] - s.loc[both, "mean_chainage_m"],
                       "se_net_change_m": np.hypot(s.loc[both, "se_chainage_m"],
                                                   e.loc[both, "se_chainage_m"])})
    return domain_net_change(tr)


# When a window's positions actually fall: count, mean date, months from the final state, share per year
def window_dates(period, name, w, model_end):
    tids = window_means(w).index
    dates = pd.concat([filter_dates(load_timeseries(str(timeseries_file(t))), *w)["date"]
                       for t in tids])
    years = dates.dt.year.value_counts(normalize=True).sort_index()
    mean = dates.dt.tz_localize(None).mean()
    return dict(period=period, end_window=name, first=w[0], last=w[1], n_positions=len(dates),
                mean_date=mean.date().isoformat(),
                months_from_model_end=round((mean - model_end).days / 30.44, 1),
                **{f"share_{y}": round(v, 3) for y, v in years.items()})


def model_change(rd, years):
    m = np.load(next(Path(rd).glob("*_shoreline_matrix.npy")))
    assert m.shape[0] == years + 1, f"{rd}: {m.shape[0] - 1} model years, expected {years}"
    return pd.Series(-(m[-1] - m[0])[D.start_real_index:D.end_real_index], index=GIS)


def score(model, target):
    m, o = smooth_like_model(model).loc[INTERIOR], target.loc[INTERIOR]
    d = m - o
    return dict(bias_m=d.mean(), rmse_m=float(np.sqrt((d ** 2).mean())),
                r=float(np.corrcoef(m, o)[0, 1]), obs_mean_m=o.mean(), model_mean_m=m.mean())


def dates_txt(w):
    a, b = map(pd.Timestamp, w)
    txt = f"{a.day} {a:%b %Y} to {b.day} {b:%b %Y}"
    return txt + "*" if b > DATA_END else txt          # * = past the data end, in the footnote


# One period's panel: the five observed ends and the run; returns the legend handles and off-axis marks
def panel(a, p, ends, dom, run, model, sc, half=POS_HALF):
    add_geographic_annotations(a, dataclasses.replace(HATTERAS_ANNOTATIONS, groin_label_y=0.03,
                                                      groin_label_side="left"))
    a.axhline(0, color=INK_MUTED, lw=0.6)
    a.set_xlim(1, 90)
    a.set_ylim(-half, half)
    a.grid(axis="y")
    open_frame(a)
    mod = smooth_like_model(model)
    off = [(f"{p['title']} model", mark_offaxis(a, GIS, mod.values, half, color="0.35"))]
    # the others first, the selected line on top
    for k in [*[k for k in ends if k != SELECTED], SELECTED]:
        _, _, c, lw = ends[k]
        y = dom[k]["net_change_lowess7_m"]
        a.plot(GIS, y, color=c, lw=lw, zorder=7 if k == SELECTED else 6)
        off.append((f"{p['title']} {k}", mark_offaxis(a, GIS, y.values, half, color=c)))
    a.plot(GIS, mod, color="0.25", lw=1.0, ls="--", zorder=8)
    a.set_ylabel("Net change (m)")
    handles = [Line2D([], [], color=c, lw=lw,
                      label=f"{name}\n{dates_txt(w)}\nRMSE {sc[k]['rmse_m']:.1f} m, "
                            f"bias {sc[k]['bias_m']:+.1f} m")
               for k, (w, name, c, lw) in ends.items()]
    handles.append(Line2D([], [], color="0.25", lw=1.0, ls="--",
                          label=f"CASCADE, {run}\nmodel state on 1 Jan {p['last_year'] + 1}\n "))
    return handles, off


def footnote(fig, ends):
    if any(pd.Timestamp(w[1]) > DATA_END for w, _, _, _ in ends.values()):
        fig.text(0.01, 0.005, f"* CoastSat ends {DATA_END.day} {DATA_END:%b %Y}; the window holds "
                 "data up to then.", fontsize=6.5, color=INK_MUTED, ha="left", va="bottom")


def figure(p, ends, dom, run, model, sc, png):
    apply_style()
    fig, a = plt.subplots(figsize=figsize("double", height=4.2), constrained_layout=True)
    handles, off = panel(a, p, ends, dom, run, model, sc)
    s, y = p["start"], p["last_year"]
    a.set_title(f"{p['title']}, {run}: model vs five end windows (start mean {s[0]} to {s[1]})",
                fontsize=9)
    a.set_xlabel(DOMAIN_AXIS_LABEL)
    fig.legend(handles=handles, loc="outside lower center", ncol=3, fontsize=7, frameon=False,
               columnspacing=2.0)
    footnote(fig, ends)
    save(fig, png, dpi=300, close=True)
    lines = " ".join(f"{name}: end mean {w[0]} to {w[1]}." for w, name, _, _ in ends.values())
    record_caption(png, (
        f"Modelled and observed net shoreline change per GIS domain, {p['title']}, {run}. "
        f"Observed: five end windows, each minus the CoastSat mean over {s[0]} to {s[1]}; "
        f"red (drawn on top) is the selected end (the last modelled year, {y}, ±6 months). {lines} CoastSat ends "
        f"{DATA_END.date()}, so windows reaching past it hold only what precedes it. Dashed: the "
        f"{run} run, {p['model_start']} to 1 Jan {y + 1}; RMSE and bias (m, interior GIS 2-89) in "
        "the legend. All at 7-domain LOWESS, GIS 1-10 raw. Shaded bands: communities; dashed lines: "
        "village centres; hatched: shoal zones; dash-dot: piers; dotted: Buxton groin. Seaward "
        "positive." + offaxis_clause(off, POS_HALF, "m")))


# Both periods on one page, calibration above test; each panel keeps its own legend (dates, scores)
def both_figure(run, per, png):
    apply_style()
    fig, axes = plt.subplots(2, 1, sharex=True, figsize=figsize("double", height=8.6),
                             constrained_layout=True)
    off = []
    for i, (a, (period, d)) in enumerate(zip(axes, per.items())):
        p = PERIODS[period]
        handles, o = panel(a, p, d["ends"], d["dom"], run, d["model"], d["sc"], half=BOTH_HALF)
        off += o
        # the same run from the ±6-month island offset, second dashed line; both scored on the selected end
        six = smooth_like_model(d["model_6mo"])
        a.plot(GIS, six, color=SIX_MO_INK, lw=1.1, ls=SIX_MO_DASH, zorder=8)
        off.append((f"{p['title']} ±6-month model", mark_offaxis(a, GIS, six.values, BOTH_HALF,
                                                                  color=SIX_MO_INK)))
        tgt = d["dom"][SELECTED]["net_change_lowess7_m"]
        s1, s6 = score(d["model"], tgt), score(d["model_6mo"], tgt)
        handles[-1].set_label(f"CASCADE, ±1 yr island offset\nstate on 1 Jan {p['last_year'] + 1}\n"
                              f"vs selected: RMSE {s1['rmse_m']:.1f} m, bias {s1['bias_m']:+.1f} m")
        handles.append(Line2D([], [], color=SIX_MO_INK, lw=1.1, ls=SIX_MO_DASH,
                              label=f"CASCADE, ±6 month island offset\nstate on 1 Jan {p['last_year'] + 1}\n"
                                    f"vs selected: RMSE {s6['rmse_m']:.1f} m, bias {s6['bias_m']:+.1f} m"))
        s = p["start"]
        a.set_title(f"({'ab'[i]}) {p['title']}: model state on 1 Jan {p['last_year'] + 1} vs five "
                    f"end windows (start mean {s[0]} to {s[1]})", fontsize=9, loc="left")
        note = (f"* CoastSat ends {DATA_END.day} {DATA_END:%b %Y}; windows marked * hold data up to then"
                if any(pd.Timestamp(w[1]) > DATA_END for w, _, _, _ in d["ends"].values()) else None)
        leg = a.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.04 if i == 0 else -0.16),
                       ncol=3, fontsize=6.5, frameon=False, columnspacing=2.0, title=note,
                       title_fontsize=6.5)
        leg.get_title().set_color(INK_MUTED)
    axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)
    fig.suptitle(f"End-window options, {run}: same rule and colour in both periods", fontsize=10)
    save(fig, png, dpi=300, close=True)
    record_caption(png, (
        f"The five end-window options in both periods, {run}. (a) Calibration 1996-2009, last "
        "modelled year 2008; (b) test 2009-2025, last modelled year 2025. Each observed line is an "
        "end mean minus the period's DEM-centred start mean, and each option uses the same rule and "
        "colour in both periods: red (drawn on top), selected, the last modelled year ±6 months; "
        "green, the last modelled year only; teal, the last two modelled years; purple, the last "
        "modelled year ±1 yr (3-yr window); orange, the model's final state ±1 yr. Long dashes: the "
        "run from the ±1-yr island offset (offset v2, the production matrix); short grey dashes: the "
        "same run from the ±6-month island offset (offset v3, 2026-10-09-6mo-starts), with the same "
        "end rates. RMSE and bias (m, interior GIS 2-89) in each legend; the two runs are scored on "
        f"the selected end. CoastSat ends {DATA_END.date()}. "
        "All at 7-domain LOWESS, GIS 1-10 raw. Seaward positive." + offaxis_clause(off, BOTH_HALF, "m")))


def main():
    pd.set_option("display.width", 250)
    both = {}                                       # run -> period -> what the two-panel figure needs
    for period, p in PERIODS.items():
        dest = HERE / period
        y = p["last_year"]
        model_end = pd.Timestamp(f"{y + 1}-01-01")
        years = y + 1 - int(p["model_start"][-4:])
        ends = ends_for(y)
        out = pd.DataFrame(index=GIS)
        dom = {}
        for k, (w, _, _, _) in ends.items():
            dom[k] = observed(p["start"], w)
            out[f"{k}_m"] = dom[k]["net_change_m"]
            out[f"{k}_lowess7_m"] = dom[k]["net_change_lowess7_m"]
        dates = pd.DataFrame([window_dates(period, k, w, model_end) for k, (w, _, _, _) in ends.items()])
        dates = dates[[c for c in dates.columns if not c.startswith("share_")]
                      + sorted(c for c in dates.columns if c.startswith("share_"))]
        dest.mkdir(parents=True, exist_ok=True)
        dates.to_csv(dest / "window_dates.csv", index=False)
        rows = []
        for run, rd in p["runs"].items():
            preset = run.split(",")[0].split(" ")[0]
            model = model_change(rd, years)
            out[f"model_{preset}_lowess7_m"] = smooth_like_model(model)
            sc = {k: score(model, dom[k]["net_change_lowess7_m"]) for k in ends}
            rows += [dict(period=period, run=run, end_window=k, **v) for k, v in sc.items()]
            figure(p, ends, dom, run, model, sc, dest / "figures" / f"end_options_{period}_{preset}.png")
            six = model_change(SIX_MO / rd.relative_to(RAW), years)
            both.setdefault(run, {})[period] = dict(ends=ends, dom=dom, model=model, sc=sc,
                                                    model_6mo=six)
        out.round(3).to_csv(dest / "domain_net_change.csv", index_label="gis_domain")
        t = pd.DataFrame(rows).round(2)
        t.to_csv(dest / "scores.csv", index=False)
        print(f"\n=== {period} ===\n" + dates.to_string(index=False) + "\n\n" + t.to_string(index=False))
        i = out.loc[INTERIOR]
        for k in list(ends)[1:]:
            d = i[f"{k}_lowess7_m"] - i[f"{SELECTED}_lowess7_m"]
            print(f"{k} minus selected, LOWESS-7 interior: mean {d.mean():+.2f}, "
                  f"RMS {np.sqrt((d ** 2).mean()):.2f}, range {d.min():+.1f}..{d.max():+.1f}")
    for run, per in both.items():
        preset = run.split(",")[0].split(" ")[0]
        both_figure(run, per, HERE / "figures" / f"end_options_both_periods_{preset}.png")


if __name__ == "__main__":
    main()
