#!/usr/bin/env python3
"""
How much does each period's target move if its end line is centred elsewhere?

    python target_end_window.py  ->  <period>/domain_net_change.csv, scores.csv, window_dates.csv,
                                     <period>/figures/end_lines_vs_model_<period>_<preset>.png
                                     for <period> = calibration, test

Calibration 1996-2009, all from the 1996 start mean (1995-10-12 to 1997-10-12):
  dem_centred  the plan's end, 2008-08-17 to 2010-08-17 (centred on the 2009 DEM)
  mid_2009     2007-07-01 to 2009-06-30 (calendar 2008 +/-6 months, no data after June 2009)
  jan_2009     2008-01-01 to 2010-01-01 (+/-1 yr of the model's end, 1 Jan 2009)
  year_2009    2009-01-01 to 2009-12-31 (every position in the end year; min 5 positions, not 10)
Test 2009-2025, all from the 2009 start mean (2008-08-17 to 2010-08-17):
  dem_2025     the plan's end, 2025-02-17 to 2026-02-17 (17 Aug 2025, the 2009 DEM date's anniversary, +/-6 months)
  cal_2024     2023-07-01 to 2025-06-30 (calendar 2024 +/-6 months)
  jan_2025     2024-01-01 to 2026-01-01 (+/-1 yr of the model's end, 1 Jan 2025)
  jan_2025_6mo 2024-07-01 to 2025-07-01 (+/-6 months of the model's end)
  year_2025    2025-01-01 to 2025-12-31 (every position in the end year)
The existing runs are scored against each, interior GIS 2-89, both sides 7-domain LOWESS.
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
RAW = REPO / "output" / "raw_runs"
GIS = pd.RangeIndex(1, 91)
INTERIOR = slice(2, 89)
POS_HALF = 100.0                                # m, the net-change figures' fixed axis
NUMBER = {2: "two", 3: "three", 4: "four", 5: "five", 6: "six"}
CAL_START, _ = NET_CHANGE_WINDOWS[(1996, 2009)]
CAL_DEM_END = ("2008-08-17", "2010-08-17")         # the target this study tested; the live one
                                                   # became calendar 2008 +/-6 months on 10-09
TEST_START, TEST_DEM_END = NET_CHANGE_WINDOWS[(2009, 2025)]
# Per period: start window, end windows (key -> dates, legend label, colour), model end, runs.
# A colour means the same kind of window in both periods, and LEGEND_SLOTS gives it the same
# legend position; labels read "centre ± half-width". Blue DEM date, red mid-year before the end, orange and purple the model end
# +/-1 yr and +/-6 months, green the end calendar year.
LEGEND_SLOTS = ["#2166ac", "#b2182b", "#e08214", "#762a83", "#1b7837", "model"]
PERIODS = {
    "calibration": dict(
        title="Calibration 1996-2009", start=CAL_START, start_key="start_1996",
        model_start="1 Jan 1996", model_end=pd.Timestamp("2009-01-01"),
        ends={
            "dem_centred": (CAL_DEM_END, "Old target: 17 Aug 2009 (DEM date) ±1 yr", "#2166ac"),
            "mid_2009": (("2007-07-01", "2009-06-30"), "Current target: 1 Jul 2008 ±1 yr", "#b2182b"),
            "jan_2009": (("2008-01-01", "2010-01-01"), "1 Jan 2009 (model end) ±1 yr", "#e08214"),
            "year_2009": (("2009-01-01", "2009-12-31"), "Calendar year 2009 (min 5 positions)", "#1b7837"),
        },
        runs={
            "zeroBE, no groin": RAW / "matrix/1996_2009/zeroBE/HAT_1996_2009_zeroBE_offsetmetres_road_bdm_nogroin",
            "edgeBE, no groin": RAW / "matrix/1996_2009/edgeBE/HAT_1996_2009_edgeBE_offsetmetres_road_bdm_nogroin",
            "domainBE + groin": RAW / "matrix/1996_2009/domainBE/HAT_1996_2009_domainBE_offsetmetres_road_bdm_groinblock",
        }),
    "test": dict(
        title="Test 2009-2025", start=TEST_START, start_key="start_2009",
        model_start="1 Jan 2009", model_end=pd.Timestamp("2025-01-01"),
        ends={
            "dem_2025": (TEST_DEM_END, "Current target: 17 Aug 2025 (DEM date) ±6 months", "#2166ac"),
            "cal_2024": (("2023-07-01", "2025-06-30"), "1 Jul 2024 ±1 yr", "#b2182b"),
            "jan_2025": (("2024-01-01", "2026-01-01"), "1 Jan 2025 (model end) ±1 yr", "#e08214"),
            "jan_2025_6mo": (("2024-07-01", "2025-07-01"), "1 Jan 2025 (model end) ±6 months", "#762a83"),
            "year_2025": (("2025-01-01", "2025-12-31"), "Calendar year 2025", "#1b7837"),
        },
        runs={
            "zeroBE, no groin": RAW / "matrix/2009_2025/zeroBE/HAT_2009_2025_zeroBE_offsetmetres_road_bdm_nourish_nogroin",
            "edgeBE, no groin": RAW / "matrix/2009_2025/edgeBE/HAT_2009_2025_edgeBE_offsetmetres_road_bdm_nourish_nogroin",
            "domainBE + groin": RAW / "matrix/2009_2025/domainBE/HAT_2009_2025_domainBE_offsetmetres_road_bdm_nourish_groinblock",
        }),
}
# -----------------------------------------------------------------------------


# Observed change between a start mean and one end window, per transect and per domain
def observed(start_w, end_w):
    s, e = window_means(start_w), window_means(end_w)
    both = s.index.intersection(e.index)
    tr = pd.DataFrame({"domain_number": s.loc[both, "domain_number"].astype(int),
                       "end_mean_chainage_m": e.loc[both, "mean_chainage_m"],
                       "net_change_m": e.loc[both, "mean_chainage_m"] - s.loc[both, "mean_chainage_m"],
                       "se_net_change_m": np.hypot(s.loc[both, "se_chainage_m"],
                                                   e.loc[both, "se_chainage_m"])})
    return tr, domain_net_change(tr)


# When a window's positions actually fall: mean date, and the share from each calendar year
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


def model_change(rd):
    m = np.load(next(Path(rd).glob("*_shoreline_matrix.npy")))
    return pd.Series(-(m[-1] - m[0])[D.start_real_index:D.end_real_index], index=GIS)


def score(model, target):
    m, o = smooth_like_model(model).loc[INTERIOR], target.loc[INTERIOR]
    d = m - o
    return dict(bias_m=d.mean(), rmse_m=float(np.sqrt((d ** 2).mean())),
                r=float(np.corrcoef(m, o)[0, 1]), obs_mean_m=o.mean(), model_mean_m=m.mean())


def figure(p, dom, run, model, png):
    apply_style()
    fig, a = plt.subplots(figsize=figsize("double", height=4.0), constrained_layout=True)
    add_geographic_annotations(a, dataclasses.replace(HATTERAS_ANNOTATIONS, groin_label_y=0.03,
                                                      groin_label_side="left"))
    a.axhline(0, color=INK_MUTED, lw=0.6)
    a.set_xlim(1, 90)
    a.set_ylim(-POS_HALF, POS_HALF)
    a.grid(axis="y")
    open_frame(a)
    mod = smooth_like_model(model)
    off = [("model", mark_offaxis(a, GIS, mod.values, POS_HALF, color="0.35"))]
    for k, (_, _, c) in p["ends"].items():
        y = dom[k]["net_change_lowess7_m"]
        a.plot(GIS, y, color=c, lw=1.4, zorder=6)
        off.append((k, mark_offaxis(a, GIS, y.values, POS_HALF, color=c)))
    a.plot(GIS, mod, color="0.35", lw=1.0, ls="--", zorder=5)
    s, n = p["start"], NUMBER[len(p["ends"])]
    a.set_title(f"{p['title']}, {run}: model vs {n} end lines (start {s[0]} to {s[1]})", fontsize=9)
    a.set_ylabel("Net change (m)")
    a.set_xlabel(DOMAIN_AXIS_LABEL)

    def dates(w):
        return " to ".join(f"{t.day} {t:%b %Y}" for t in map(pd.Timestamp, w))
    end_txt = f"{p['model_end'].day} {p['model_end']:%b %Y}"
    handles = {c: Line2D([], [], color=c, lw=1.4, label=f"{lab}\n({dates(w)})")
               for w, lab, c in p["ends"].values()}
    handles["model"] = Line2D([], [], color="0.35", lw=1.0, ls="--",
                              label=f"CASCADE, {run}\nmodel state on {end_txt}")
    blank = Line2D([], [], lw=0, label=" \n ")             # keeps an unused slot's place
    fig.legend(handles=[handles.get(k, blank) for k in LEGEND_SLOTS],
               loc="outside lower center", ncol=3, fontsize=7, frameon=False, columnspacing=2.5)
    save(fig, png, dpi=300, close=True)
    lines = " ".join(f"{lab}: end mean {w[0]} to {w[1]}." for w, lab, _ in p["ends"].values())
    record_caption(png, (
        f"Modelled and observed {p['title'].split()[-1]} net shoreline change per GIS domain, "
        f"{p['title'].split()[0].lower()} period, {run}. Observed: {n} end lines, each minus the "
        f"CoastSat mean over {s[0]} to {s[1]}. {lines} Dashed: the {run} run, {p['model_start']} "
        f"to {end_txt}. All at 7-domain LOWESS, GIS 1-10 raw. Shaded bands: "
        "communities; dashed lines: village centres; hatched: shoal zones; dash-dot: piers; "
        "dotted: Buxton groin. Seaward positive." + offaxis_clause(off, POS_HALF, "m")))


def main():
    pd.set_option("display.width", 250)
    for period, p in PERIODS.items():
        dest = HERE / period
        out = pd.DataFrame(index=GIS)
        rows, summary = [], []
        ends = {k: w for k, (w, _, _) in p["ends"].items()}
        current = next(iter(ends))
        trs, dom = {}, {}
        for k, w in ends.items():
            trs[k], dom[k] = observed(p["start"], w)
            out[f"{k}_m"] = dom[k]["net_change_m"]
            out[f"{k}_lowess7_m"] = dom[k]["net_change_lowess7_m"]
        date_rows = [window_dates(period, k, w, p["model_end"]) for k, w in ends.items()]
        date_rows.append(window_dates(period, p["start_key"], p["start"], p["model_end"]))

        for k in list(ends)[1:]:
            both = trs[k].index.intersection(trs[current].index)
            d = (trs[k].loc[both, "end_mean_chainage_m"]
                 - trs[current].loc[both, "end_mean_chainage_m"])
            out[f"diff_{k}_m"] = out[f"{k}_m"] - out[f"{current}_m"]
            out[f"diff_{k}_lowess7_m"] = out[f"{k}_lowess7_m"] - out[f"{current}_lowess7_m"]
            summary.append((k, current, d))

        for run, rd in p["runs"].items():
            preset = run.split(",")[0].split(" ")[0]
            model = model_change(rd)
            out[f"model_{preset}_lowess7_m"] = smooth_like_model(model)
            for k in ends:
                rows.append(dict(period=period, run=run, end_window=k,
                                 **score(model, dom[k]["net_change_lowess7_m"])))
            figure(p, dom, run, model, dest / "figures" / f"end_lines_vs_model_{period}_{preset}.png")

        out.round(3).to_csv(dest / "domain_net_change.csv", index_label="gis_domain")
        dates = pd.DataFrame(date_rows).dropna(axis=1, how="all")
        dates = dates[[c for c in dates.columns if not c.startswith("share_")]
                      + sorted(c for c in dates.columns if c.startswith("share_"))]
        dates.to_csv(dest / "window_dates.csv", index=False)
        t = pd.DataFrame(rows).round(2)
        t.to_csv(dest / "scores.csv", index=False)

        print(f"\n=== {period} ===\n" + dates.to_string(index=False))
        i = out.loc[INTERIOR]
        for k, cur, d in summary:
            print(f"\n{k} end minus {cur} end, per transect ({len(d)}): mean {d.mean():+.2f}, "
                  f"median {d.median():+.2f}, SD {d.std():.2f}, |d|>10 m: {(d.abs() > 10).sum()}")
            r = i[f"diff_{k}_m"]
            print(f"  per domain interior: raw diff mean {r.mean():+.2f}, RMS "
                  f"{np.sqrt((r ** 2).mean()):.2f}, range {r.min():+.1f}..{r.max():+.1f}, "
                  f"seaward in {(r > 0).mean():.0%}")
            l = i[f"diff_{k}_lowess7_m"]
            print(f"  LOWESS-7 interior: diff mean {l.mean():+.2f}, RMS {np.sqrt((l ** 2).mean()):.2f}, "
                  f"range {l.min():+.1f}..{l.max():+.1f}, seaward in {(l > 0).mean():.0%}; r between targets "
                  f"{np.corrcoef(i[f'{cur}_lowess7_m'], i[f'{k}_lowess7_m'])[0, 1]:.3f}")
            print(f"  ends: GIS 1 diff {out.loc[1, f'diff_{k}_m']:+.1f}, "
                  f"GIS 90 diff {out.loc[90, f'diff_{k}_m']:+.1f}")
        print("\n" + t.to_string(index=False))


if __name__ == "__main__":
    main()
