"""
dsas_vs_coastsat_datematched.py
==============================================================================
DSAS against CoastSat with the CoastSat side anchored on the SHORELINE SURVEY
DATES rather than on a calendar window (Hannah, 2026-09-22, "do both": at
+/-30 days and at +/-6 months).

    4-comparisons/dsas_vs_coastsat/survey_dates/dsas_vs_coastsat_datematched.png
    4-comparisons/dsas_vs_coastsat/survey_dates/slides/ (--slide)
    4-comparisons/dsas_vs_coastsat/survey_dates/supporting/*.csv

THE METHOD
    Not an OLS over a window. For each CoastSat transect, the mean shoreline
    position within +/-W days of the START survey date is subtracted from the
    mean within +/-W days of the END date and divided by the interval -- the
    same endpoint method coastsat_endpoint.py uses against the dune line, run
    through the same endpoint_by_transect(). That is what "match the imagery
    dates" means: both sources then describe motion between the same two
    moments, not between two calendar years.

THE ARCHIVE DID NOT DO THIS
    5-scr/archive/coastsat_lrr_superseded_20260810/dsas_coastsat_specific_dates/
    is named for this method but was produced with SURVEY_DATES = [] in
    coastsat_domain_lrr_specific_dates.py, so it fell through to the
    continuous-range mode. Its median n_obs (417) matches the plain calendar
    fit (414.5) and its statistics are the calendar comparison's to two
    decimals. This script is the analysis that folder's name promised.

THE TWO ANCHORS, AND ONE DISCREPANCY
    1997-09-27 and 2019-09-07. The 2019 date is what nc_shorelines.geojson
    carries in SHR_DATE (Wet-Dry, 32 features). The 1997 date is Hannah's
    (2026-09-22) and matches the date commented into
    coastsat_domain_lrr_specific_dates.py, so two independent records of the
    survey agree on it.

    THE INVENTORY AGREED ONLY AFTER IT WAS FIXED. nc_shorelines.geojson
    stamped every 1997 feature 1/1/1997, a placeholder Hannah entered when the
    date was not to hand. On 2026-09-22 the two features covering this study
    area -- "Outer Banks - National Seashore" and "Outer Banks - North of
    Oregon Inlet" -- were restamped 9/27/1997; the other 21, elsewhere in the
    state and flown on other days, still carry the placeholder. See
    1-observations/shoreline_inventory/PROVENANCE.md.

    THAT EDIT DOES NOT REACH THIS SCRIPT. The dates below are module
    constants; nothing here opens the geojson, so the restamp changed no
    number in this comparison (verified by re-running: identical to three
    decimals). The two records now agree, which is worth having, but they are
    still two records.

    Both window widths are still drawn, because the two anchors are late
    September and early September -- nearly the same point in the seasonal
    cycle, so a tight window is meaningful, and agreement between the widths
    says the result does not depend on how much of the year is swept in.

    The two DSAS ends are different proxies (MHW in 1997, wet/dry in 2019)
    while CoastSat is one proxy throughout. That is a property of the DSAS
    rate, not of the matching, and no window width fixes it.

USAGE
    python scripts/input_prep/5-scr/4-comparisons/dsas_vs_coastsat/dsas_vs_coastsat_datematched.py
    python scripts/input_prep/5-scr/4-comparisons/dsas_vs_coastsat/dsas_vs_coastsat_datematched.py --slide
==============================================================================
"""

from __future__ import annotations

import argparse
import datetime as dt
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_REPO = next(_p for _p in Path(__file__).resolve().parents
             if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(_REPO / "scripts"))
sys.path.insert(0, str(_REPO / "scripts" / "input_prep" / "5-scr" / "lib"))
import scr_paths  # noqa: E402,F401  (5-scr sibling modules onto sys.path)

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import MultipleLocator  # noqa: E402

from coastsat_vs_duneline import (DAYS_PER_YEAR, endpoint_by_transect,  # noqa: E402
                                  load_chainage)
from site_layer import hat_observed_rates as obs  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C_1984, C_1997, DOMAIN_AXIS_LABEL, INK_MUTED, apply_style, caption,
    figsize, open_frame, save, structures, support_dir, town_bands,
)

N = 90
# The survey date, from Hannah (2026-09-22), matching the date commented into
# coastsat_domain_lrr_specific_dates.py. The inventory's 1997-01-01 is her own
# placeholder, not an alternative reading; see the header.
START_DATE = dt.datetime(1997, 9, 27, tzinfo=dt.timezone.utc)
END_DATE = dt.datetime(2019, 9, 7, tzinfo=dt.timezone.utc)
HALF_DAYS = (30.0, 182.0)          # "do both"
DSAS_WINDOW = (1997, 2019)
OUT = obs.COMPARISONS / "dsas_vs_coastsat" / "survey_dates"
STEM = "dsas_vs_coastsat_datematched"
C_DSAS, C_CS = C_1984, C_1997
LW, LW_SLIDE = 1.1, 0.8
SLIDE_W_IN, SLIDE_H_IN = 3.4, 4.2


def dsas():
    df = pd.read_csv(obs.DSAS_ROOT / f"dsas_{DSAS_WINDOW[0]}_{DSAS_WINDOW[1]}_rates.csv")
    df = df.rename(columns={"domain_id": "domain", "MEAN_LRR": "lrr"})
    return (df[["domain", "lrr"]].groupby("domain")["lrr"].mean()
            .reindex(range(1, N + 1)))


def coastsat_datematched(half_days, lookup, cache):
    """Per-domain endpoint rate from CoastSat, anchored on the survey dates.

    Returns (per-domain Series, per-transect frame). A transect with no
    position inside one of the two windows yields NaN and drops out.
    """
    ep = endpoint_by_transect(lookup, START_DATE, END_DATE, half_days, cache)
    per_domain = (ep.dropna(subset=["endpoint_rate_m_yr"])
                  .groupby("domain_number")["endpoint_rate_m_yr"].mean()
                  .reindex(range(1, N + 1)))
    return per_domain, ep


def agreement(a, b):
    ok = a.notna() & b.notna()
    d = (b[ok] - a[ok]).to_numpy(float)
    r = (float(np.corrcoef(a[ok], b[ok])[0, 1]) if ok.sum() > 2 else float("nan"))
    return dict(n=int(ok.sum()), bias=float(d.mean()),
                rmse=float(np.sqrt((d ** 2).mean())), r=r)


def figure(d, cs, stats, cover, slide=False):
    x = np.arange(1, N + 1)
    lw = LW_SLIDE if slide else LW
    size = (figsize(SLIDE_W_IN, height=SLIDE_H_IN) if slide
            else figsize("double", height=5.0))
    fig, axes = plt.subplots(len(HALF_DAYS), 1, sharex=True, sharey=True,
                             figsize=size, constrained_layout=True)
    for i, (w, ax) in enumerate(zip(HALF_DAYS, axes)):
        ax.set_xlim(0.5, N + 0.5)
        town_bands(ax, label=(i == 0 and not slide))
        ax.axhline(0, color=INK_MUTED, lw=0.6, zorder=2)
        ax.plot(x, d, color=C_DSAS, lw=lw, zorder=4)
        ax.plot(x, cs[w], color=C_CS, lw=lw, zorder=5)
        structures(ax, i == 0 and not slide, 4.0)
        ax.yaxis.set_major_locator(MultipleLocator(2))
        ax.xaxis.set_major_locator(MultipleLocator(20 if slide else 10))
        ax.xaxis.set_minor_locator(MultipleLocator(10 if slide else 5))
        ax.yaxis.grid(True, zorder=0)
        ax.set_axisbelow(True)
        ax.set_title(f"({chr(97 + i)}) ±{_w(w)}", loc="left", fontsize=8)
        open_frame(ax)
    axes[-1].set_xlabel("GIS domain (S → N)" if slide else DOMAIN_AXIS_LABEL)
    fig.supylabel("Shoreline change rate (m/yr)", fontsize=9)
    fig.legend([Line2D([], [], color=C_DSAS, lw=lw),
                Line2D([], [], color=C_CS, lw=lw)],
               (["DSAS", "CoastSat"] if slide else
                ["DSAS (digitized shorelines)",
                 "CoastSat (endpoint at the survey dates)"]),
               loc="outside lower center", ncol=2, frameon=False)

    years = (END_DATE - START_DATE).days / DAYS_PER_YEAR
    t = "; ".join(
        f"±{_w(w)} bias {stats[w]['bias']:+.2f}, RMSE {stats[w]['rmse']:.2f} m/yr, "
        f"r {stats[w]['r']:.2f} (n {stats[w]['n']} domains, "
        f"{cover[w]['transects']} of {cover[w]['total']} transects, median "
        f"{cover[w]['median_start']:.0f} and {cover[w]['median_end']:.0f} "
        "positions per end)" for w in HALF_DAYS)
    caption(fig, (
        "Observed shoreline change rate by GIS domain (1 at Cape Point, 90 at "
        "Pea Island), 1997–2019, with the CoastSat side anchored on the "
        "shoreline survey dates instead of a calendar window. For each CoastSat "
        "transect the mean position within ±W days of 2019-09-07 minus the mean "
        f"within ±W days of 1997-09-27, over {years:.2f} years, averaged per "
        "500 m domain; DSAS is its per-domain mean LRR over the same two "
        "shorelines. Seaward positive, no smoothing on either side. The two "
        "anchors sit within three weeks of each other in the seasonal cycle, "
        "so the tight window is meaningful; the wide one is drawn beside it to "
        "show the result does not depend on how much of the year is swept in. "
        "Note the shoreline inventory stamps its 1997 features 1997-01-01, a "
        "year stamp; the survey date used here is the one recorded with the "
        "analysis. "
        f"Agreement, CoastSat minus DSAS: {t}. The DSAS ends are also different "
        "proxies — MHW in 1997, wet/dry in 2019 — while CoastSat is one proxy "
        "throughout; no window width addresses that. Village spans are shaded; "
        "the solid hairline is the Buxton groin and the dotted hairlines are "
        "the Avon and Rodanthe piers."))
    return fig


def _w(half_days):
    return "30 days" if half_days < 100 else "6 months"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--slide", action="store_true",
                    help=f"{SLIDE_W_IN:g} in canvas, to <stem>_slide.png")
    slide = ap.parse_args(argv).slide
    apply_style()

    lookup = pd.read_csv(obs.transect_lookup())
    lookup = lookup[lookup["domain_number"].between(1, N)]
    cache: dict = {}
    d = dsas()
    cs, stats, cover = {}, {}, {}
    cols = {"domain": np.arange(1, N + 1), "dsas_1997_2019": d.to_numpy(float)}
    for w in HALF_DAYS:
        series, ep = coastsat_datematched(w, lookup, cache)
        cs[w] = series
        stats[w] = agreement(d, series)
        cover[w] = dict(total=len(ep),
                        transects=int(ep["endpoint_rate_m_yr"].notna().sum()),
                        median_start=float(ep["n_start"].median()),
                        median_end=float(ep["n_end"].median()))
        tag = "30d" if w < 100 else "6mo"
        cols[f"coastsat_{tag}"] = series.to_numpy(float)
        cols[f"difference_{tag}"] = (series - d).to_numpy(float)
        v, c = stats[w], cover[w]
        print(f"±{_w(w):9s}: n {v['n']:2d} domains  bias {v['bias']:+.2f}  "
              f"RMSE {v['rmse']:.2f}  r {v['r']:.2f}  | transects with both "
              f"ends {c['transects']}/{c['total']}, median positions per end "
              f"{c['median_start']:.0f}/{c['median_end']:.0f}")

    fig = figure(d, cs, stats, cover, slide=slide)
    folder = OUT / "slides" if slide else OUT
    out = save(fig, folder / (STEM + ("_slide" if slide else "")))
    plt.close(fig)
    table = support_dir(OUT) / f"{STEM}.csv"
    pd.DataFrame(cols).round(3).to_csv(table, index=False)
    for p in out + [table]:
        print(f"wrote    {p.relative_to(_REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
