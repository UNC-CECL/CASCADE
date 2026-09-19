"""
coastsat_lrr_projected.py
==============================================================================
The long-term shoreline change rate projected to a DISTANCE, beside the
distance the shoreline actually moved. Built 2026-09-19 (Hannah, by
interview, for her advisor's "total change in shoreline position from the
long-term rate").

PROJECTED   per CoastSat transect, its 1996-2024 LRR (3-rates/coastsat/lrr,
            the OLS slope through every satellite position in the calendar
            window) times the window length, end - start = 28 yr: where the
            shoreline would be if it had moved at its long-term rate the
            whole time. Per domain the mean of its transects, which is the
            domain mean LRR times 28, so the alongshore PATTERN is the LRR
            figure's; only the unit changes.

OBSERVED    per transect, the mean position over the whole END calendar year
            minus the mean over the whole START calendar year (all of 2024
            minus all of 1996). Both means are centred mid-year, so the span
            is the same 28 yr as the projection, and both use only data inside
            the LRR's window. Not the dune-date endpoint in
            3-rates/coastsat/endpoint (1997-10 to 2023-07, 25.7 yr), which is
            a shorter span.

    observed - projected is how far the actual change departs from the
    long-term trend: positive where the shoreline ended up more seaward than
    its trend predicts (a fill, an accelerating accretion), negative where
    more landward. SEAWARD IS POSITIVE throughout, as in every 3-rates product.

OUTPUT   data/hatteras_init/5-scr/3-rates/coastsat/lrr_projected/<start>_<end>/
    transect_lrr_projected.csv         per transect: lrr_m_yr, its uncertainty,
                                       projected_change_m (and its uncertainty),
                                       the two calendar-year means with their
                                       counts, observed_change_m,
                                       observed_minus_projected_m
    domain_lrr_projected_summary.csv   per domain: the means of those, std,
                                       pct_landward of each
    lrr_projected_<start>_<end>.png    projected as the house-style fill and
                                       dots, observed as a black line; PDF and
                                       caption under supporting/
    PROVENANCE.md

USAGE
    python scripts/input_prep/5-scr/coastsat_lrr_projected/coastsat_lrr_projected.py
    python ... --windows 1996_2024 2010_2024
==============================================================================
"""

from __future__ import annotations

import argparse
import datetime as dt
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_REPO = next(_p for _p in Path(__file__).resolve().parents
             if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(_REPO / "scripts"))
for _sub in ("duneline_vs_coastsat", "", "CoastSat"):
    sys.path.insert(0, str(_REPO / "scripts" / "input_prep" / "5-scr" / _sub))

from duneline_vs_coastsat import load_chainage  # noqa: E402
import rates_figures as rf  # noqa: E402  (the 3-rates drawing helpers)
from rates_figures import cw, plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.legend_handler import HandlerTuple  # noqa: E402
from site_layer.hat_figure_style import INK, apply_style, caption, save  # noqa: E402
from site_layer.hat_observed_rates import (  # noqa: E402
    COASTSAT_LRR_PROJECTED_ROOT, COASTSAT_LRR_ROOT,
)

N_DOMAINS = 90
WINDOWS = [(1996, 2024)]
TRANSECT_FILE = "transect_lrr_projected.csv"
DOMAIN_FILE = "domain_lrr_projected_summary.csv"
OBS_LW = 1.1


def year_mean(df, year):
    """(mean position, n) over one calendar year, or (nan, 0)."""
    sel = df.loc[df["date"].dt.year == year, "chainage"]
    return (float(sel.mean()), int(sel.size)) if sel.size else (np.nan, 0)


def build(start: int, end: int, cache: dict) -> dict:
    years = end - start
    lrr = pd.read_csv(COASTSAT_LRR_ROOT / f"{start}_{end}" / "transect_lrr_full.csv")
    lrr = lrr[lrr["domain_number"].between(1, N_DOMAINS)].copy()
    lrr["domain_number"] = lrr["domain_number"].astype(int)

    rows = []
    for tid in lrr["transect_id"]:
        if tid not in cache:
            cache[tid] = load_chainage(tid)
        df = cache[tid]
        m0, n0 = year_mean(df, start) if df is not None else (np.nan, 0)
        m1, n1 = year_mean(df, end) if df is not None else (np.nan, 0)
        rows.append(dict(transect_id=tid, position_start_m=m0, n_start=n0,
                         position_end_m=m1, n_end=n1))
    t = lrr[["transect_id", "domain_number", "lrr_m_yr", "unc_m_yr", "r_squared",
             "n_obs", "start_date", "end_date"]].merge(pd.DataFrame(rows), on="transect_id")
    t["projected_change_m"] = t["lrr_m_yr"] * years
    t["projected_unc_m"] = t["unc_m_yr"] * years
    t["observed_change_m"] = t["position_end_m"] - t["position_start_m"]
    t["observed_minus_projected_m"] = t["observed_change_m"] - t["projected_change_m"]
    t = t.assign(window=f"{start}_{end}", projection_years=years)

    g = t.groupby("domain_number")
    both = t[t["observed_change_m"].notna() & t["projected_change_m"].notna()].groupby("domain_number")
    dom = pd.DataFrame({
        "n_transects": g.size(),
        "mean_lrr_m_yr": g["lrr_m_yr"].mean(),
        "mean_projected_change_m": g["projected_change_m"].mean(),
        "std_projected_change_m": g["projected_change_m"].std(),
        "pct_projected_landward": g["projected_change_m"].apply(lambda s: 100.0 * (s < 0).mean()),
        "n_observed": both.size(),
        "mean_observed_change_m": both["observed_change_m"].mean(),
        "std_observed_change_m": both["observed_change_m"].std(),
        "pct_observed_landward": both["observed_change_m"].apply(lambda s: 100.0 * (s < 0).mean()),
        "mean_observed_minus_projected_m": both["observed_minus_projected_m"].mean(),
    }).round(3).reindex(range(1, N_DOMAINS + 1)).rename_axis("domain_number").reset_index()

    # The domain projection must be the LRR table's own domain mean x years.
    ref = pd.read_csv(COASTSAT_LRR_ROOT / f"{start}_{end}" / "domain_lrr_summary.csv")
    ref = ref.set_index(ref["domain_number"].astype(int))["mean_lrr"]
    check = float(np.nanmax(np.abs(dom.set_index("domain_number")["mean_lrr_m_yr"]
                                   - ref.reindex(range(1, N_DOMAINS + 1)))))

    out = COASTSAT_LRR_PROJECTED_ROOT / f"{start}_{end}"
    out.mkdir(parents=True, exist_ok=True)
    t.to_csv(out / TRANSECT_FILE, index=False, float_format="%.4f")
    dom.to_csv(out / DOMAIN_FILE, index=False)
    return dict(start=start, end=end, years=years, t=t, dom=dom, out=out, check=check)


def figure(r) -> list:
    s, e, years, t, dom = r["start"], r["end"], r["years"], r["t"], r["dom"]
    ext = np.nanmax(np.abs(np.r_[dom["mean_projected_change_m"], dom["mean_observed_change_m"]]))
    half = float(math.ceil((ext + 5) / rf.Y_STEP_M) * rf.Y_STEP_M)
    tick = 10.0 if half <= 60 else 20.0
    tt, x = rf._along(t)
    fig, ax, n_out = rf._draw(rf._frame(dom, "mean_projected_change_m"), x,
                              tt["projected_change_m"].to_numpy(float), half,
                              "Net change in shoreline position (m)",
                              tick, cw.fills_in(s, e), std=False)
    ax.plot(dom["domain_number"], dom["mean_observed_change_m"], color=INK,
            lw=OBS_LW, zorder=12, solid_capstyle="round")   # over the pier labels' boxes
    h = [(Line2D([], [], color=cw.C_ACCRETE, lw=1.0), Line2D([], [], color=cw.C_ERODE, lw=1.0)),
         (Line2D([], [], color=cw.C_ACCRETE, marker="o", ms=2.2, lw=0),
          Line2D([], [], color=cw.C_ERODE, marker="o", ms=2.2, lw=0)),
         Line2D([], [], color=INK, lw=OBS_LW)]
    labels = [f"Projected change (LRR × {years} yr)",
              "Projected change, individual transects",
              "Observed change (endpoint)"]
    fig.legend(h, labels, loc="outside lower center", ncol=2, frameon=False,
               handler_map={tuple: HandlerTuple(ndivide=None, pad=0.3)})
    ok = dom.dropna(subset=["mean_observed_change_m", "mean_projected_change_m"])
    diff = ok["mean_observed_minus_projected_m"]
    caption(fig, (
        f"Projected net change in shoreline position by GIS domain (1 at Cape Point, "
        f"90 at Pea Island), {s}–{e}: for each CoastSat transect the {s}–{e} linear "
        "regression rate (the ordinary-least-squares slope through every satellite "
        f"position from 1 January {s} to 31 December {e}) multiplied by {years} yr, "
        "the distance the shoreline would have moved at its long-term rate. The "
        "coloured line and fill are the mean of the ~10 transects in each 500 m "
        "domain, blue seaward and red landward, so the alongshore pattern is the "
        f"{s}–{e} rate's; the dots are the single transects, blue or red by their own "
        "sign"
        + (f" ({n_out} beyond the axis, drawn as open circles at its edge)" if n_out else "")
        + f". The black line is the OBSERVED ENDPOINT change over the same {years} yr: per "
        f"transect the mean position over all of {e} minus the mean over all of {s}, "
        "averaged per domain. Where it lies above the fill the shoreline ended more "
        "seaward than its trend predicts; below, more landward. Over the "
        f"{len(ok)} domains with both, observed minus projected averages "
        f"{diff.mean():+.1f} m (range {diff.min():+.1f} to {diff.max():+.1f} m). "
        "Seaward positive, in metres. " + rf._marks_clause(s, e)
        + f" The y axis is ±{half:g} m."))
    out = save(fig, r["out"] / f"lrr_projected_{s}_{e}")
    plt.close(fig)
    return out


def provenance(r) -> None:
    s, e, years, t, dom = r["start"], r["end"], r["years"], r["t"], r["dom"]
    ok = dom.dropna(subset=["mean_observed_change_m"])
    diff = ok["mean_observed_minus_projected_m"]
    corr = np.corrcoef(ok["mean_projected_change_m"], ok["mean_observed_change_m"])[0, 1]
    no_obs = int(t["observed_change_m"].isna().sum())
    (r["out"] / "PROVENANCE.md").write_text("\n".join([
        f"# 3-rates/coastsat/lrr_projected/{s}_{e} - provenance",
        "",
        f"Written {dt.datetime.now():%Y-%m-%d %H:%M} by "
        "scripts/input_prep/5-scr/coastsat_lrr_projected/coastsat_lrr_projected.py.",
        "",
        f"**Projected** = the transect's {s}-{e} LRR (`../../lrr/{s}_{e}/transect_lrr_full.csv`) "
        f"x {years} yr ({e} - {s}). Per domain the mean over its transects; it equals the "
        f"LRR table's `mean_lrr` x {years} to {r['check']:.1e} m/yr.",
        "",
        f"**Observed** = mean CoastSat position over calendar {e} minus the mean over "
        f"calendar {s}, per transect, from the same time series the LRR is fitted to. "
        f"Median positions per transect: {t['n_start'].median():.0f} in {s}, "
        f"{t['n_end'].median():.0f} in {e}. {no_obs} of {len(t)} transects have no "
        "position in one of the two years and no observed change.",
        "",
        "Seaward positive. `observed_minus_projected_m` > 0 means the shoreline ended "
        "more seaward than its long-term trend predicts.",
        "",
        "## Island summary",
        "",
        f"Domain mean projected {dom['mean_projected_change_m'].mean():+.1f} m, observed "
        f"{ok['mean_observed_change_m'].mean():+.1f} m. Projected landward in "
        f"{int((dom['mean_projected_change_m'] < 0).sum())} of {N_DOMAINS} domains, observed "
        f"in {int((ok['mean_observed_change_m'] < 0).sum())} of {len(ok)}. Observed minus "
        f"projected per domain: mean {diff.mean():+.1f} m, range {diff.min():+.1f} to "
        f"{diff.max():+.1f} m; r(projected, observed) = {corr:.2f}.",
        "",
    ]), encoding="utf-8")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="LRR projected to a distance, beside the observed change.")
    ap.add_argument("--windows", nargs="+", metavar="START_END")
    a = ap.parse_args(argv)
    wins = ([tuple(int(x) for x in w.split("_")) for w in a.windows] if a.windows else WINDOWS)
    apply_style()
    cache: dict = {}
    for s, e in wins:
        r = build(s, e, cache)
        figure(r)
        provenance(r)
        d = r["dom"]
        print(f"{s}_{e}  x{r['years']} yr  projected {d['mean_projected_change_m'].mean():+6.1f} m  "
              f"observed {d['mean_observed_change_m'].mean():+6.1f} m  "
              f"obs-proj {d['mean_observed_minus_projected_m'].mean():+5.1f} m  "
              f"(check {r['check']:.1e})  -> {r['out'].relative_to(_REPO)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
