#!/usr/bin/env python3
"""
Two calendar-year ends for the test target, beside the DEM-anniversary end in use.

    python test_end_options.py  ->  domain_net_change.csv, scores.csv, window_dates.csv,
                                    figures/test_end_options_edgeBE.png

Test 2009-2025, all from the 2009 start mean (2008-08-17 to 2010-08-17, DEM-centred, kept):
  dem_2025   in use: 2025-02-17 to 2026-02-17 (17 Aug 2025, the 2009 DEM date's anniversary, ±6 months)
  cal_2025   calendar 2025 ±6 months: 2024-07-01 to 2026-06-30 (calibration's rule; CoastSat ends 2026-01-13)
  year_2025  calendar year 2025 only: 2025-01-01 to 2025-12-31
The model runs 2009 through 2025 (17 yr, final state 1 Jan 2026). The existing matrix runs are
scored against each, interior GIS 2-89, both sides 7-domain LOWESS. No new runs.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-10
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = next(p for p in HERE.parents if (p / "pyproject.toml").exists())
NC = REPO / "scripts" / "input_prep" / "5-scr" / "3-rates" / "coastsat" / "net_change"
MS = REPO / "scripts" / "input_prep" / "5-scr" / "1-observations" / "mean_shoreline"
sys.path[:0] = [str(REPO / "scripts"), str(NC), str(MS)]
from coastsat_net_change import domain_net_change, smooth_like_model, window_means  # noqa: E402
from coastsat_mean_shoreline import timeseries_file  # noqa: E402
from coastsat_lrr import filter_dates, load_timeseries  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    DOMAIN_AXIS_LABEL, INK_MUTED, _title, apply_style, figsize, mark_offaxis,
    offaxis_clause, record_caption, save, town_bands)
from site_layer.hat_observed_rates import NET_CHANGE_WINDOWS  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_DOMAINS as D  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
RAW = REPO / "output" / "raw_runs" / "matrix" / "2009_2025"
GIS = pd.RangeIndex(1, 91)
INTERIOR = slice(2, 89)
HALF = 60.0                                     # m, axis half-range
START, DEM_END = NET_CHANGE_WINDOWS[(2009, 2025)]
MODEL_END = pd.Timestamp("2026-01-01")          # final state of the 17-yr test run
DATA_END = "13 Jan 2026"                        # last CoastSat position
# key -> (dates, legend label, colour); blue and green as in the end-window figures
ENDS = {
    "dem_2025": (DEM_END, "In use: 17 Aug 2025 (DEM date) ±6 months", "#2166ac"),
    "cal_2025": (("2024-07-01", "2026-06-30"), "Option 1: calendar 2025 ±6 months", "#35978f"),
    "year_2025": (("2025-01-01", "2025-12-31"), "Option 2: calendar year 2025 only", "#1b7837"),
}
RUNS = {
    "zeroBE, no groin": RAW / "zeroBE/HAT_2009_2025_zeroBE_offsetmetres_road_bdm_nourish_nogroin",
    "edgeBE, no groin": RAW / "edgeBE/HAT_2009_2025_edgeBE_offsetmetres_road_bdm_nourish_nogroin",
    "domainBE + groin": RAW / "domainBE/HAT_2009_2025_domainBE_offsetmetres_road_bdm_nourish_groinblock",
}
FIG_RUN = "edgeBE, no groin"
PNG = HERE / "figures" / "test_end_options_edgeBE.png"
# -----------------------------------------------------------------------------


def observed(start_w, end_w):
    s, e = window_means(start_w), window_means(end_w)
    both = s.index.intersection(e.index)
    tr = pd.DataFrame({"domain_number": s.loc[both, "domain_number"].astype(int),
                       "net_change_m": e.loc[both, "mean_chainage_m"] - s.loc[both, "mean_chainage_m"],
                       "se_net_change_m": np.hypot(s.loc[both, "se_chainage_m"],
                                                   e.loc[both, "se_chainage_m"])})
    return domain_net_change(tr)["net_change_lowess7_m"]


# When a window's positions actually fall: count, mean date, months from the model's final state
def window_dates(name, w):
    tids = window_means(w).index
    dates = pd.concat([filter_dates(load_timeseries(str(timeseries_file(t))), *w)["date"]
                       for t in tids])
    mean = dates.dt.tz_localize(None).mean()
    return dict(end_window=name, first=w[0], last=w[1], n_positions=len(dates),
                last_position=dates.max().date().isoformat(), mean_date=mean.date().isoformat(),
                months_from_model_end=round((mean - MODEL_END).days / 30.44, 1))


def model_change(rd):
    m = np.load(next(Path(rd).glob("*_shoreline_matrix.npy")))
    assert m.shape[0] == 18, f"{rd} is not the 17-yr test run"
    return pd.Series(-(m[-1] - m[0])[D.start_real_index:D.end_real_index], index=GIS)


def score(model, target):
    m, o = smooth_like_model(model).loc[INTERIOR], target.loc[INTERIOR]
    d = m - o
    return dict(bias_m=round(d.mean(), 2), rmse_m=round(float(np.sqrt((d ** 2).mean())), 2),
                r=round(float(np.corrcoef(m, o)[0, 1]), 2), obs_mean_m=round(o.mean(), 2),
                model_mean_m=round(m.mean(), 2))


def dates_txt(w):
    return " to ".join(f"{t.day} {t:%b %Y}" for t in map(pd.Timestamp, w))


def figure(tgt, mod, sc):
    x = GIS.to_numpy()
    apply_style()
    fig, axes = plt.subplots(len(ENDS) + 1, 1, sharex=True, sharey=True,
                             figsize=figsize("double", height=8.5), constrained_layout=True)
    off = []

    def frame(a):
        a.axhline(0, color=INK_MUTED, lw=0.6)
        a.set_ylim(-HALF, HALF)
        a.set_xlim(x.min(), x.max())
        a.set_xticks(range(5, 91, 5), minor=True)
        a.grid(axis="y")
        a.grid(axis="x", which="both", color="0.93", lw=0.5)
        town_bands(a, where="top", strip=0.07, label=True)
        a.set_ylabel("Net change (m)")

    a = axes[0]
    frame(a)
    for k, (w, lab, c) in ENDS.items():
        a.plot(x, tgt[k], color=c, lw=1.4, label=f"{lab} ({dates_txt(w)})")
        off.append((f"{lab} target", mark_offaxis(a, x, tgt[k].values, HALF, color=c)))
    a.legend(loc="lower center", bbox_to_anchor=(0.58, 0.0), fontsize=6.5, frameon=False)
    _title(a, 0, "The three test targets, each from the 2009 start mean (17 Aug 2008 to 17 Aug 2010)")

    for i, (k, (w, lab, c)) in enumerate(ENDS.items(), start=1):
        a = axes[i]
        frame(a)
        s = sc[k]
        a.plot(x, tgt[k], color=c, lw=1.4, label=f"Observed, end {dates_txt(w)}")
        a.plot(x, mod, color="0.25", lw=1.4, ls="--", label=f"CASCADE {FIG_RUN}, state on 1 Jan 2026")
        leg = a.legend(loc="lower center", bbox_to_anchor=(0.58, 0.0), fontsize=6.5, frameon=False,
                       title=f"bias {s['bias_m']:+.1f} m, RMSE {s['rmse_m']:.1f} m, r {s['r']:.2f}",
                       title_fontsize=6.5)
        leg.get_title().set_color(INK_MUTED)
        _title(a, i, lab)
    off.append((f"{FIG_RUN} run", mark_offaxis(axes[1], x, mod.values, HALF, color="0.25")))
    axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)
    save(fig, PNG, dpi=300, close=True)
    record_caption(PNG, (
        "Three candidate end windows for the test target (2009-2025), each minus the CoastSat mean "
        "over the 2009 start window (17 Aug 2008 to 17 Aug 2010). (a) The three observed targets: "
        "blue, the end in use, 17 Aug 2025 ±6 months; teal, calendar 2025 ±6 months (1 Jul 2024 to "
        f"30 Jun 2026, calibration's rule; CoastSat ends {DATA_END}, so the window is off-centre); "
        "green, calendar year 2025 only. (b-d) Each target against the edgeBE run (no groin), whose "
        "final state is 1 Jan 2026; scores (m, interior GIS 2-89) in each legend. The edgeBE end "
        "rates were solved on the end in use, so GIS 1 and 90 favour (b). All at 7-domain LOWESS. "
        "Seaward positive." + offaxis_clause(off, HALF, "m")))
    print(f"wrote {PNG}")


def main():
    pd.set_option("display.width", 200)
    tgt = pd.DataFrame({k: observed(START, w) for k, (w, _, _) in ENDS.items()}, index=GIS)
    out = tgt.add_suffix("_lowess7_m")
    rows, fig_scores, fig_model = [], {}, None
    for run, rd in RUNS.items():
        model = model_change(rd)
        preset = run.split(",")[0].split(" ")[0]
        out[f"model_{preset}_lowess7_m"] = smooth_like_model(model)
        for k in ENDS:
            s = score(model, tgt[k])
            rows.append(dict(run=run, end_window=k, **s))
            if run == FIG_RUN:
                fig_scores[k] = s
        if run == FIG_RUN:
            fig_model = smooth_like_model(model)
    out.round(3).to_csv(HERE / "domain_net_change.csv", index_label="gis_domain")
    sc = pd.DataFrame(rows)
    sc.to_csv(HERE / "scores.csv", index=False)
    wd = pd.DataFrame([window_dates(k, w) for k, (w, _, _) in ENDS.items()])
    wd.to_csv(HERE / "window_dates.csv", index=False)
    print(wd.to_string(index=False))
    print("\n" + sc.to_string(index=False))
    i = tgt.loc[INTERIOR]
    for k in list(ENDS)[1:]:
        d = i[k] - i["dem_2025"]
        print(f"{k} minus dem_2025, interior LOWESS-7: mean {d.mean():+.2f}, RMS "
              f"{np.sqrt((d ** 2).mean()):.2f}, range {d.min():+.1f}..{d.max():+.1f}")
    figure(tgt, fig_model, fig_scores)


if __name__ == "__main__":
    main()
