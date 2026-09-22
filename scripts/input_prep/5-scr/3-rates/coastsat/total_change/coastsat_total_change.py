"""
coastsat_total_change.py
==============================================================================
The CoastSat linear regression rate turned into a DISTANCE, beside the
distance the shoreline actually moved. Built 2026-09-19 (Hannah, by
interview, for her advisor's "total change in shoreline position from the
long-term rate"); split into two named products 2026-09-21 (Hannah, by
interview) after the folder called `lrr_projected/` turned out to hold no
projections at all.

THE VOCABULARY.  A rate turned into a distance is named by the window it was
FITTED on, never by the arithmetic:

  TOTAL SHORELINE CHANGE   --product total  ->  3-rates/coastsat/total_change/
      The rate is evaluated over the SAME window it was fitted on.
      LRR(1996-2010) x 14 yr, LRR(2010-2024) x 14 yr, LRR(1996-2024) x 28 yr.
      Nothing is extrapolated, so nothing is projected. This is what the
      whole of the old `lrr_projected/` tree actually was.

  PROJECTED SHORELINE      --product projected  ->  3-rates/coastsat/projected/
  CHANGE
      The 1996-2024 rate carried onto a window it was NOT fitted on:
      LRR(1996-2024) x 14 yr over 1996-2010, and the same over 2010-2024.
      1996_2024 is deliberately absent -- there it would BE the total change.
      This is the pairing the model's CoastSat target uses in both halves
      (output/comparisons/target_comparison/projected/), here on the
      observations alone.

  OBSERVED CHANGE          both products, unchanged
      No rate anywhere: per transect, the mean position over the whole END
      calendar year minus the mean over the whole START calendar year (all of
      2010 minus all of 1996). Both means are centred mid-year, so the span
      is the same as the rate's multiply, and both use only data inside the
      window. Not the dune-date endpoint in 3-rates/coastsat/endpoint
      (1997-10 to 2023-07, 25.7 yr), which is a shorter span.

    observed - (total or projected) is how far the actual change departs from
    the trend: positive where the shoreline ended up more seaward than the
    trend predicts, negative where more landward. Under --product projected
    it is the more interesting residual of the two, because the rate there
    was never fitted to the window it is being judged over. SEAWARD IS
    POSITIVE throughout, as in every 3-rates product.

SMOOTHED    the same comparison after an alongshore LOESS (Hannah, by
            interview, 2026-09-21). The rate the MODEL is graded against is not
            the raw rate: it is raw over GIS 1-10 and a 10-domain LOESS of the
            transect values beyond (cascade_pipeline.coastsat_loess). The raw
            comparison above therefore tests the fairness of a quantity nobody
            uses; this one tests the target as it is actually applied.

            Note LOESS commutes with the x years multiply -- the weights depend
            only on the transect positions and the robust reweighting is scale
            equivariant -- so smoothing the RATE and smoothing the DISTANCE
            give the same number to machine precision. Nothing here turns on
            the order; what matters is that BOTH sides are smoothed, at the
            same window, so the residual is not a smoothed quantity minus an
            unsmoothed one.

            Windows 3, 5 and 10 domains (1.5, 2.5, 5.0 km) are all built. The
            sweep is the point: if the residual collapses as the window
            widens, the departures from trend are transect-scale estimation
            noise; if a departure survives 10 domains, the rate genuinely
            fails there, at the scale the model resolves. Read the bias and
            the RMS residual, NOT r -- smoothing strips high-frequency
            variance that is uncorrelated between the two sides, so r rises
            whether or not the smoothing is telling the truth.

FIGURE TITLES carry quantity, window and method, so a figure pulled out of
its folder still says which of the two products it is (Hannah, 2026-09-21):
"Total shoreline change, 1996-2010 (CoastSat LRR 1996-2010 x 14 yr)" against
"Projected shoreline change, 1996-2010 (CoastSat LRR 1996-2024 x 14 yr)" --
the method names the window the rate was FITTED on, so the reader never has
to trust the folder.

OUTPUT   data/hatteras_init/5-scr/3-rates/coastsat/<product>/<start>_<end>/
    transect_<product>.csv             per transect: lrr_m_yr, its uncertainty,
                                       <product>_change_m (and its uncertainty),
                                       the two calendar-year means with their
                                       counts, observed_change_m,
                                       observed_minus_<total|projected>_m
    domain_<product>_summary.csv       per domain: the means of those, std,
                                       pct_landward of each
    <product>_<start>_<end>.png        the rate's distance as the house-style
                                       fill and dots, observed as a black
                                       line; PDF and caption under supporting/
    PROVENANCE.md
    smoothed/<product>_smoothed_<start>_<end>_w<NN>.png
                                       one per LOESS window, same axis as each
                                       other so the windows can be read side by
                                       side; PDFs and captions under
                                       smoothed/supporting/
    smoothed/tables/domain_smoothed.csv        long: one row per domain per
                                       window (0 = the raw product above)
    smoothed/tables/residual_by_scale.csv      the sweep: bias, RMS, range,
                                       sign agreement and r per window
    smoothed/PROVENANCE.md

USAGE
    python scripts/input_prep/5-scr/3-rates/coastsat/total_change/coastsat_total_change.py
    python ... --product projected
    python ... --product both
    python ... --windows 1996_2024 2010_2024
    python ... --smooth-windows 3 5 10 | --no-smoothed
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
sys.path.insert(0, str(_REPO / "scripts" / "input_prep" / "5-scr" / "lib"))
import scr_paths  # noqa: E402,F401  (5-scr sibling modules onto sys.path)

from coastsat_vs_duneline import load_chainage  # noqa: E402
import rates_figures as rf  # noqa: E402  (the 3-rates drawing helpers)
from rates_figures import cw, plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.legend_handler import HandlerTuple  # noqa: E402
from matplotlib.ticker import MultipleLocator  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C, DOMAIN_AXIS_LABEL, INK, SMOOTH_RAMP, apply_style, caption,
    compare_header, figsize, mark_offaxis, offaxis_clause, save,
)
from site_layer.hat_observed_rates import (  # noqa: E402
    COASTSAT_LRR_ROOT, COASTSAT_PROJECTED_ROOT, COASTSAT_TOTAL_CHANGE_ROOT,
    PROJECTED_RATE_WINDOW, WINDOW_ROLE,
)
# The model target's own smoother, imported rather than re-implemented so the
# 10-domain figure here IS the treatment the runs are graded under.
from cascade_pipeline.coastsat_loess import spliced_loess_series  # noqa: E402
from cascade_pipeline.domains import DEFAULT_DOMAINS  # noqa: E402

N_DOMAINS = 90
# The canonical 1996 -> 2010 -> 2024 chain: the full period and its two halves.
# Which of these a product is defined for is on the Product below, because
# `projected` has no 1996_2024 (see PROJECTED).
OBS_LW = 1.1
# The observed CoastSat line, PURPLE since 2026-09-22 (Hannah). It was
# black, and a black line in 4-comparisons/shoreline_vs_duneline is the
# DUNE LINE -- same glyph, two meanings across trees, which is exactly
# how this one got read as the dune line. The house ACCENT purple, so it
# is a colour the project already uses rather than a new one.
C_OBSERVED = C["ACCENT"]
# ONE fixed metre axis on every figure here (Hannah, 2026-09-22), shared
# with 4-comparisons/shoreline_vs_duneline/total_change and
# output/comparisons/target_comparison, so a figure from any of the three
# can be laid beside another without rescaling by eye. It is FIXED, not a
# floor: a window whose data exceeds it is marked at the edge and named in
# the caption (hat_figure_style.mark_offaxis) rather than given its own
# axis, which would defeat the point. Only Cape Point does, in practice.
Y_HALF_M = 100.0
Y_TICK_M = 20.0
NL = chr(10)          # the provenance writers join on it


class Product:
    """One of the two named products. The ONLY thing that differs between them
    is which window the rate is read from; everything downstream -- the
    observed side, the figures, the LOESS sweep -- is identical, which is the
    point of building both from one script.

    Attributes:
        key: folder name under 3-rates/coastsat/ and the file stem.
        noun: how the quantity is named in a title, caption or legend.
        tok: the token that replaces `rate` in the written column names, so
            every CSV says which product it is without its path.
        windows: the change windows this product is defined for.
    """

    def __init__(self, key, noun, tok, root, windows, rate_window, method):
        self.key, self.noun, self.tok = key, noun, tok
        self.root, self.windows = root, windows
        self._rate_window, self._method = rate_window, method

    def rate_window(self, s, e):
        """The window the LRR is FITTED on, which is what names the product."""
        return self._rate_window(s, e)

    def method(self, s, e):
        """The parenthetical in a figure title: WHERE THE RATE CAME FROM and
        how the distance was made -- "CoastSat LRR 1996-2010 x 14 yr".

        The fit window is in the string on purpose (Hannah, 2026-09-21): the
        reader compares it against the window in the title, and the two being
        equal or not IS the difference between total change and a projection.
        So the two products read in parallel and differ in one number:
            Total shoreline change, 1996-2010 (CoastSat LRR 1996-2010 x 14 yr)
            Projected shoreline change, 1996-2010 (CoastSat LRR 1996-2024 x 14 yr)
        """
        return self._method(s, e)

    def cols(self, df):
        """Internal `rate` column names -> the product's written ones."""
        keep = {"mean_lrr_m_yr", "lrr_window"}
        ren = {c: (f"{self.tok}_change_m" if c == "rate_m" else c.replace("rate", self.tok))
               for c in df.columns if "rate" in c and c not in keep}
        return df.rename(columns=ren)

    @property
    def transect_file(self):
        return f"transect_{self.key}.csv"

    @property
    def domain_file(self):
        return f"domain_{self.key}_summary.csv"


# Each window's own rate over its own years. Nothing is extrapolated.
TOTAL = Product(
    "total_change", "Total shoreline change", "total", COASTSAT_TOTAL_CHANGE_ROOT,
    [(1996, 2024), (1996, 2010), (2010, 2024)],
    rate_window=lambda s, e: (s, e),
    method=lambda s, e: f"CoastSat LRR {s}–{e} × {e - s} yr",
)
# The 1996-2024 rate carried onto a window it was not fitted on. 1996_2024 is
# absent on purpose: there the rate window IS the change window, so the answer
# is TOTAL, and building it here would put the same numbers under two names.
PROJECTED = Product(
    "projected", "Projected shoreline change", "projected", COASTSAT_PROJECTED_ROOT,
    [(1996, 2010), (2010, 2024)],
    rate_window=lambda s, e: PROJECTED_RATE_WINDOW,
    method=lambda s, e: (f"CoastSat LRR {PROJECTED_RATE_WINDOW[0]}–"
                         f"{PROJECTED_RATE_WINDOW[1]} × {e - s} yr"),
)
PRODUCTS = {p.key: p for p in (TOTAL, PROJECTED)}

# LOESS window widths in domain units (1 domain = 500 m). 10 is the model
# target's window (coastsat_loess.LoessConfig.window_domains); 3 and 5 are
# there to show how fast the residual collapses with scale.
SMOOTH_WINDOWS = (3, 5, 10)
# GIS 1..SPLICE_DOMAINS keep their raw domain means instead of the LOESS --
# coastsat_loess.LoessConfig.skip_southern_domains, the boundary treatment at
# Oregon Inlet. Applied to the OBSERVED side too, so the two never differ in
# treatment at any domain.
SPLICE_DOMAINS = 10


def year_mean(df, year):
    """(mean position, n) over one calendar year, or (nan, 0)."""
    sel = df.loc[df["date"].dt.year == year, "chainage"]
    return (float(sel.mean()), int(sel.size)) if sel.size else (np.nan, 0)


def build(start: int, end: int, cache: dict, prod: Product = TOTAL) -> dict:
    """The product's rate turned into a distance over `start`-`end`, beside the
    observed change over the same years.

    The rate is read from `prod.rate_window(start, end)`, which is the window
    it was FITTED on -- the same window for TOTAL, always 1996-2024 for
    PROJECTED. That single line is the whole difference between the two
    products; everything below is shared.
    """
    years = end - start
    rs, re_ = prod.rate_window(start, end)
    lrr = pd.read_csv(COASTSAT_LRR_ROOT / f"{rs}_{re_}" / "transect_lrr_full.csv")
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
    # Internal names are neutral (`rate`); Product.cols() renames them to the
    # product's own on the way out, so no CSV is ambiguous about which it is.
    t["rate_change_m"] = t["lrr_m_yr"] * years
    t["rate_unc_m"] = t["unc_m_yr"] * years
    t["observed_change_m"] = t["position_end_m"] - t["position_start_m"]
    t["observed_minus_rate_m"] = t["observed_change_m"] - t["rate_change_m"]
    t = t.assign(window=f"{start}_{end}", lrr_window=f"{rs}_{re_}", span_years=years)

    g = t.groupby("domain_number")
    both = t[t["observed_change_m"].notna() & t["rate_change_m"].notna()].groupby("domain_number")
    dom = pd.DataFrame({
        "n_transects": g.size(),
        "mean_lrr_m_yr": g["lrr_m_yr"].mean(),
        "mean_rate_change_m": g["rate_change_m"].mean(),
        "std_rate_change_m": g["rate_change_m"].std(),
        "pct_rate_landward": g["rate_change_m"].apply(lambda s: 100.0 * (s < 0).mean()),
        "n_observed": both.size(),
        "mean_observed_change_m": both["observed_change_m"].mean(),
        "std_observed_change_m": both["observed_change_m"].std(),
        "pct_observed_landward": both["observed_change_m"].apply(lambda s: 100.0 * (s < 0).mean()),
        "mean_observed_minus_rate_m": both["observed_minus_rate_m"].mean(),
    }).round(3).reindex(range(1, N_DOMAINS + 1)).rename_axis("domain_number").reset_index()

    # The domain distance must be the LRR table's own domain mean x years.
    ref = pd.read_csv(COASTSAT_LRR_ROOT / f"{rs}_{re_}" / "domain_lrr_summary.csv")
    ref = ref.set_index(ref["domain_number"].astype(int))["mean_lrr"]
    check = float(np.nanmax(np.abs(dom.set_index("domain_number")["mean_lrr_m_yr"]
                                   - ref.reindex(range(1, N_DOMAINS + 1)))))

    out = prod.root / f"{start}_{end}"
    out.mkdir(parents=True, exist_ok=True)
    prod.cols(t).to_csv(out / prod.transect_file, index=False, float_format="%.4f")
    prod.cols(dom).to_csv(out / prod.domain_file, index=False)
    return dict(start=start, end=end, years=years, t=t, dom=dom, out=out, check=check,
                prod=prod, rate_window=(rs, re_))


def figure(r) -> list:
    s, e, years, t, dom = r["start"], r["end"], r["years"], r["t"], r["dom"]
    prod, (rs, re_) = r["prod"], r["rate_window"]
    half, tick = Y_HALF_M, Y_TICK_M
    tt, x = rf._along(t)
    fig, ax, n_out = rf._draw(rf._frame(dom, "mean_rate_change_m"), x,
                              tt["rate_change_m"].to_numpy(float), half,
                              "Net change in shoreline position (m)",
                              tick, cw.fills_in(s, e), std=False)
    ax.plot(dom["domain_number"], dom["mean_observed_change_m"], color=C_OBSERVED,
            lw=OBS_LW, zorder=12, solid_capstyle="round")   # over the pier labels' boxes
    # Both series are clipped by the fixed axis, so say what left it.
    off = [(f"the {prod.tok} change",
            mark_offaxis(ax, dom["domain_number"], dom["mean_rate_change_m"],
                         half, color=INK)),
           ("the observed change",
            mark_offaxis(ax, dom["domain_number"], dom["mean_observed_change_m"],
                         half, color=C_OBSERVED))]
    h = [(Line2D([], [], color=cw.C_ACCRETE, lw=1.0), Line2D([], [], color=cw.C_ERODE, lw=1.0)),
         (Line2D([], [], color=cw.C_ACCRETE, marker="o", ms=2.2, lw=0),
          Line2D([], [], color=cw.C_ERODE, marker="o", ms=2.2, lw=0)),
         Line2D([], [], color=C_OBSERVED, lw=OBS_LW)]
    labels = [f"{prod.noun} ({prod.method(s, e)})",
              f"{prod.noun}, individual transects",
              "CoastSat observed change (calendar-year endpoints)"]
    fig.legend(h, labels, loc="outside lower center", ncol=2, frameon=False,
               handler_map={tuple: HandlerTuple(ndivide=None, pad=0.3)})
    # Quantity, window, method (Hannah, 2026-09-21). The method is what tells
    # total from projected at a glance -- "LRR × 14 yr" against "1996–2024 LRR
    # × 14 yr" -- so it is on the canvas, not left to the folder. draw_fills
    # puts its bars at 1.025 in axes fractions with the year above them, so
    # the title has to clear those when the window contains a fill.
    ax.set_title(f"{prod.noun}, {s}–{e} ({prod.method(s, e)})",
                 pad=20 if cw.fills_in(s, e) else 6)
    # Both series here are CoastSat. 4-comparisons is where a CoastSat
    # series meets a dune-line one; this tree never mixes sources, and
    # after the observed line was read as the dune line it says so; it is
        # purple now for the same reason.
    compare_header(fig, [
        f"{s}–{e}   ·   BOTH series are CoastSat — no dune line on this figure",
        f"{prod.noun.lower()}: {prod.method(s, e)}   ·   observed: CoastSat mean position, all of {e} minus all of {s}"])
    ok = dom.dropna(subset=["mean_observed_change_m", "mean_rate_change_m"])
    diff = ok["mean_observed_minus_rate_m"]
    caption(fig, (
        f"{prod.noun} by GIS domain (1 at Cape Point, 90 at Pea Island), {s}–{e}: "
        f"for each CoastSat transect the {rs}–{re_} linear regression rate (the "
        "ordinary-least-squares slope through every satellite position from "
        f"1 January {rs} to 31 December {re_}) multiplied by {years} yr, the distance "
        "the shoreline would have moved at that rate. "
        + (f"The rate is fitted on {rs}–{re_} and evaluated over the SAME window, so "
           "nothing is extrapolated — this is total change, not a projection. "
           if (rs, re_) == (s, e) else
           f"The rate is fitted on {rs}–{re_} but evaluated over {s}–{e}, a window it "
           "was not fitted on, which is what makes this a PROJECTION: it is the "
           f"long-term trend asked what {s}–{e} should have looked like. ")
        + "The coloured line and fill are the mean of the ~10 transects in each 500 m "
        "domain, blue seaward and red landward, so the alongshore pattern is the "
        f"{rs}–{re_} rate's; the dots are the single transects, blue or red by their "
        "own sign"
        + (f" ({n_out} beyond the axis, drawn as open circles at its edge)" if n_out else "")
        + f". The PURPLE line is the OBSERVED ENDPOINT change over the same {years} yr: per "
        f"transect the mean position over all of {e} minus the mean over all of {s}, "
        "averaged per domain — no rate anywhere in it. Where it lies above the fill "
        "the shoreline ended more seaward than the trend implies; below, more "
        f"landward. Over the {len(ok)} domains with both, observed minus "
        f"{prod.tok} averages {diff.mean():+.1f} m (range {diff.min():+.1f} to "
        f"{diff.max():+.1f} m). "
        "Seaward positive, in metres. " + rf._marks_clause(s, e)
        + f" The y axis is ±{half:g} m, fixed, the same on every metre figure"
        " in 3-rates, 4-comparisons and target_comparison."
        + offaxis_clause(off, half)))
    out = save(fig, r["out"] / f"{prod.key}_{s}_{e}")
    plt.close(fig)
    return out


def _smooth_series(dom_ids, along_m, values, window):
    """One alongshore LOESS pass at transect resolution, averaged to domains,
    with GIS 1..SPLICE_DOMAINS put back to their raw domain means -- the
    scoring target's own two steps, shared with the smoothing-scale sweep in
    analyze_output/compare_runs/HAT_smoothing_scale.py.

    Args:
        dom_ids, along_m, values: per-transect domain id, along-coast distance
            in metres, and the quantity to smooth (projected or observed).
        window: LOESS window width in domain units.

    Returns:
        (Series indexed 1..N_DOMAINS, the lowess frac used).
    """
    return spliced_loess_series(dom_ids, along_m, values, window,
                                skip=SPLICE_DOMAINS)


def smooth(r, windows=SMOOTH_WINDOWS) -> dict:
    """The rate-vs-observed comparison repeated under the model target's
    alongshore LOESS, at each window in `windows`. BOTH sides get the same
    pass and the same splice, so no window compares a smoothed quantity with
    an unsmoothed one. Window 0 in the output is the raw product."""
    t = r["t"].sort_values(["domain_number", "transect_id"]).reset_index(drop=True)
    tt, x_dom = rf._along(t)                                 # x in domain units, for the dots
    along_m = x_dom * DEFAULT_DOMAINS.domain_spacing_m       # lowess is given metres
    ids = tt["domain_number"].to_numpy(int)
    proj = tt["rate_change_m"].to_numpy(float)
    obs = tt["observed_change_m"].to_numpy(float)

    idx = pd.RangeIndex(1, N_DOMAINS + 1, name="domain_number")
    raw = r["dom"].set_index("domain_number")
    series = {0: (raw["mean_rate_change_m"].reindex(idx),
                  raw["mean_observed_change_m"].reindex(idx), float("nan"))}
    for w in windows:
        p, frac = _smooth_series(ids, along_m, proj, w)
        o, _ = _smooth_series(ids, along_m, obs, w)
        series[w] = (p, o, frac)

    km = DEFAULT_DOMAINS.domain_spacing_m / 1000.0
    rows, stats = [], []
    for w in [0] + list(windows):
        p, o, frac = series[w]
        pv, ov = p.to_numpy(float), o.to_numpy(float)
        rows.append(pd.DataFrame({
            "domain_number": np.asarray(idx), "window_domains": w,
            "window_km": np.nan if w == 0 else w * km,
            "rate_m": pv, "observed_m": ov, "observed_minus_rate_m": ov - pv}))
        ok = np.isfinite(pv) & np.isfinite(ov)
        a, b = pv[ok], ov[ok]
        d = b - a
        stats.append(dict(
            window_domains=w, window_km=np.nan if w == 0 else w * km,
            loess_frac=frac, n_domains=int(ok.sum()),
            bias_m=float(d.mean()), rms_residual_m=float(np.sqrt((d ** 2).mean())),
            residual_min_m=float(d.min()), residual_max_m=float(d.max()),
            sd_rate_m=float(a.std(ddof=1)), sd_observed_m=float(b.std(ddof=1)),
            pct_sign_agreement=float(100.0 * (np.sign(a) == np.sign(b)).mean()),
            # Named for what it is: a symmetric smoother strips variance that is
            # uncorrelated between the two sides, so this climbs with the window
            # whether or not the smoothing is right. It is not a score.
            r_inflated_by_smoothing=float(np.corrcoef(a, b)[0, 1])))

    out = r["out"] / "smoothed"
    (out / "tables").mkdir(parents=True, exist_ok=True)
    long = pd.concat(rows, ignore_index=True).round(3)
    st = pd.DataFrame(stats).round(3)
    prod = r["prod"]
    prod.cols(long).to_csv(out / "tables" / "domain_smoothed.csv", index=False)
    prod.cols(st).to_csv(out / "tables" / "residual_by_scale.csv", index=False)
    return dict(series=series, long=long, stats=st, out=out, windows=tuple(windows),
                x_dom=x_dom, proj=proj)


def _smooth_bounds(sm):
    """The fixed metre axis, as everywhere else here (Hannah, 2026-09-22).
    Kept as a function so the call sites read the same as before."""
    return Y_HALF_M, Y_TICK_M


def overlay_bounds(sms):
    """One y half-range and tick for EVERY window's overlay, from the
    rate-derived series alone.

    The per-window panels bound on projected AND observed together; the
    overlay draws no observed side, so bounding it that way would set the axis
    from a series that is not on the figure. These three are meant to be read
    against each other -- 28 yr against two 14 yr halves -- so they share one
    bound, and it is not the panels' bound. Said in each caption.
    """
    return Y_HALF_M, Y_TICK_M


def smooth_overlay_figure(r, sm, half=None, tick=None) -> list:
    """Every LOESS width's distance on ONE panel, no observed side (Hannah,
    2026-09-21).

    The per-window panels above each answer "does the trend hold HERE"; this
    one answers "what does the window do to the target", which needs the
    curves on top of each other and nothing else competing for the eye.

    Note this is the RATE figure of input_prep/5-scr/coastsat_lrr_smoothing_windows.py
    in metres: LOESS commutes with the x years multiply, so the curves have
    the same shape and only the units differ. It is drawn because metres is
    the unit the model and the dune line are read in, not because it shows a
    different field."""
    s, e, years = r["start"], r["end"], r["years"]
    prod, (rs, re_) = r["prod"], r["rate_window"]
    if half is None:
        half, tick = overlay_bounds([sm])
    windows = (0,) + tuple(sm["windows"])
    km_of = DEFAULT_DOMAINS.domain_spacing_m / 1000.0
    role = WINDOW_ROLE.get((s, e))

    fig, ax = plt.subplots(figsize=figsize("double", aspect=0.40),
                           constrained_layout=True)
    # mean all-NaN draws the frame, grid, village bands and structures with no
    # sign fill: four curves share the panel, so the blue/red pair is not
    # available and the ordered ramp below carries the width instead.
    frame = pd.DataFrame({"domain_number": np.arange(1, N_DOMAINS + 1),
                          "mean_lrr": np.nan, "std_lrr": 0.0})
    cw.draw_panel(ax, frame, half, std=False)
    cw.draw_shoals(ax, label=True)
    fills = cw.fills_in(s, e)
    if fills:
        cw.draw_fills(ax, fills, half)
    for i, w in enumerate(windows):
        p = sm["series"][w][0]
        x, y = np.asarray(p.index, dtype=float), p.to_numpy(float)
        colour = SMOOTH_RAMP[i % len(SMOOTH_RAMP)]
        if w:
            ax.plot(x, y, color=colour, lw=1.5, zorder=10 + i, solid_capstyle="round")
        else:
            ax.plot(x, y, color=colour, lw=0.8, marker="o", ms=1.8, zorder=10,
                    solid_capstyle="round")
    ax.yaxis.set_major_locator(MultipleLocator(tick))
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    ax.set_ylabel("Net change in shoreline position (m)")
    off = [(f"the {w * km_of:g} km curve" if w else "the unsmoothed domain means",
            mark_offaxis(ax, np.asarray(sm["series"][w][0].index, dtype=float),
                         sm["series"][w][0].to_numpy(float), half, color=INK))
           for w in windows]
    # Quantity, window, method, like every other figure in the tree since
    # 2026-09-21. draw_fills puts its bars at 1.025 in axes fractions and the
    # year above them, so the title has to clear that when the window contains
    # a fill -- at the default pad it lands on the 2022 labels. The window's
    # ROLE in the 1996-2010-2024 chain used to be the title; it is in the
    # caption now, because the product is the thing a reader cannot recover.
    ax.set_title(f"{prod.noun}, {s}–{e} ({prod.method(s, e)}, every LOESS width)",
                 pad=20 if fills else 6)

    h = [Line2D([], [], color=SMOOTH_RAMP[0], lw=0.8, marker="o", ms=2.2)]
    h += [Line2D([], [], color=SMOOTH_RAMP[i % len(SMOOTH_RAMP)], lw=1.5)
          for i, w in enumerate(windows) if w]
    labels = ["Unsmoothed domain means"] + [
        f"LOESS {w * km_of:g} km ({w} domains)" for w in windows if w]
    fig.legend(h, labels, loc="outside lower center", ncol=len(h), frameon=False)

    # How far apart the windows are, in the unit of the axis.
    spread = {w: float((sm["series"][w][0] - sm["series"][0][0]).abs().max())
              for w in windows if w}
    spread_txt = "; ".join(f"{w * km_of:g} km up to {v:.0f} m" for w, v in spread.items())
    # The rate figure exists only where coastsat_lrr_smoothing_windows.py has been run;
    # a cross-reference to a file that is not there is worse than none.
    rate_png = (COASTSAT_LRR_ROOT / f"{rs}_{re_}" / f"smoothing_windows_{rs}_{re_}.png")
    rate_ref = (
        f"This is the rate figure `3-rates/coastsat/lrr/{rs}_{re_}/{rate_png.name}` in "
        f"metres — a LOESS commutes with the × {years} yr multiply, so the curves have "
        "the same shape and only the units differ. " if rate_png.exists() else "")
    role_txt = (f"This window is the {role.lower()} of the 1996–2010–2024 chain. "
                if role else "")
    caption(fig, (
        role_txt
        + f"{prod.noun} by GIS domain (1 at Cape Point, "
        f"90 at Pea Island), {s}–{e}, at each alongshore smoothing width. Each "
        f"transect's {rs}–{re_} linear regression rate × {years} yr is the distance the "
        "shoreline would have moved at that rate"
        + ("" if (rs, re_) == (s, e) else
           f", the rate being fitted on {rs}–{re_} and carried onto {s}–{e}")
        + "; the palest line with "
        "markers is the mean of those over the ~10 transects in each 500 m domain, "
        "unsmoothed, and the three heavier curves are the same quantity after the "
        f"LOESS the model target is built through, at {', '.join(f'{w * km_of:g} km' for w in windows[1:-1])} "
        f"and {windows[-1] * km_of:g} km, light to dark, the darkest being the "
        f"{windows[-1]}-domain window every run is graded at. Every curve keeps the "
        f"raw domain means over GIS 1–{SPLICE_DOMAINS} — the boundary treatment at "
        "Oregon Inlet — so all four are identical there by construction. The observed "
        "endpoint change is deliberately not drawn: this figure is about what the "
        "window does to the target, and the observed comparison is the three panels "
        f"beside it. Smoothing moves the projection by {spread_txt} at the most "
        "affected domain, against a raw domain-mean range of "
        f"{sm['series'][0][0].min():+.0f} to {sm['series'][0][0].max():+.0f} m. "
        + rate_ref
        + "Seaward positive, in metres. {} The y axis is ±{:g} m, fixed, so the 28 yr "
        "figure reads against the two 14 yr ones and against every other metre figure. "
        "It is the fixed metre axis every figure in this tree uses."
    ).format(rf._marks_clause(s, e), half) + offaxis_clause(off, half))
    written = save(fig, sm["out"] / f"{prod.key}_smoothed_{s}_{e}_overlay")
    plt.close(fig)
    return written


def smooth_figures(r, sm) -> list:
    """One panel per window, all on the same y axis so they read side by side."""
    s, e, years = r["start"], r["end"], r["years"]
    prod, (rs, re_) = r["prod"], r["rate_window"]
    half, tick = _smooth_bounds(sm)
    by_w = sm["stats"].set_index("window_domains")
    raw_row = by_w.loc[0]
    written = []
    for w in sm["windows"]:
        p, o, frac = sm["series"][w]
        row = by_w.loc[w]
        km = w * DEFAULT_DOMAINS.domain_spacing_m / 1000.0
        dom_w = pd.DataFrame({"domain_number": np.asarray(p.index, dtype=int),
                              "rate_m": p.to_numpy(float)})
        fig, ax, n_out = rf._draw(rf._frame(dom_w, "rate_m"), sm["x_dom"], sm["proj"],
                                  half, "Net change in shoreline position (m)",
                                  tick, cw.fills_in(s, e), std=False)
        ax.plot(np.asarray(o.index, dtype=int), o.to_numpy(float), color=C_OBSERVED,
                lw=OBS_LW, zorder=12, solid_capstyle="round")
        off = [(f"the smoothed {prod.tok} change",
                mark_offaxis(ax, np.asarray(p.index, dtype=int),
                             p.to_numpy(float), half, color=INK)),
               ("the smoothed observed change",
                mark_offaxis(ax, np.asarray(o.index, dtype=int),
                             o.to_numpy(float), half, color=C_OBSERVED))]
        h = [(Line2D([], [], color=cw.C_ACCRETE, lw=1.0), Line2D([], [], color=cw.C_ERODE, lw=1.0)),
             (Line2D([], [], color=cw.C_ACCRETE, marker="o", ms=2.2, lw=0),
              Line2D([], [], color=cw.C_ERODE, marker="o", ms=2.2, lw=0)),
             Line2D([], [], color=C_OBSERVED, lw=OBS_LW)]
        labels = [f"{prod.noun}, LOESS {km:g} km ({prod.method(s, e)})",
                  f"{prod.noun}, individual transects (raw)",
                  f"CoastSat observed change, LOESS {km:g} km"]
        fig.legend(h, labels, loc="outside lower center", ncol=2, frameon=False,
                   handler_map={tuple: HandlerTuple(ndivide=None, pad=0.3)})
        # Quantity, window, method -- plus the smoothing width, which is the
        # only thing that separates these panels from each other.
        ax.set_title(f"{prod.noun}, {s}–{e} ({prod.method(s, e)}, LOESS {km:g} km)",
                     pad=20 if cw.fills_in(s, e) else 6)
        # Both series here are CoastSat. 4-comparisons is where a CoastSat
        # series meets a dune-line one; this tree never mixes sources, and
        # after the observed line was read as the dune line it says so; it is
        # purple now for the same reason.
        compare_header(fig, [
            f"{s}–{e}   ·   BOTH series are CoastSat — no dune line on this figure",
            f"{prod.noun.lower()}: {prod.method(s, e)}   ·   observed: CoastSat mean position, all of {e} minus all of {s}"])
        caption(fig, (
            f"{prod.noun} and observed change, {s}–{e}, both "
            f"passed through the same alongshore LOESS of {w} domains ({km:g} km). "
            "The quantity the model is graded against is not the raw rate but this "
            "one — raw domain means over GIS 1–10 and a LOESS of the transect values "
            f"beyond — so this panel tests the {rs}–{re_} trend as it is actually "
            f"applied, not as it is estimated. {prod.noun.upper()} (coloured line and "
            f"fill, blue seaward and red landward) is each transect's {rs}–{re_} linear "
            f"regression rate × {years} yr, smoothed alongshore at transect resolution "
            "and then averaged per 500 m domain"
            + ("" if (rs, re_) == (s, e) else
               f" — the rate is fitted on {rs}–{re_} and carried onto {s}–{e}, a window "
               "it was not fitted on")
            + f". OBSERVED (purple line) is the mean CoastSat "
            f"position over all of {e} minus the mean over all of {s}, per transect, "
            "through the identical pass — both sides are smoothed, so the gap between "
            "them is not an artefact of treating them differently. The dots are the "
            "RAW individual transects, unsmoothed, drawn so the width of the cloud the "
            "line was drawn through stays visible"
            + (f" ({n_out} beyond the axis, drawn as open circles at its edge)" if n_out else "")
            + f". Over {int(row['n_domains'])} domains the residual (observed minus "
            f"{prod.tok}) has mean {row['bias_m']:+.1f} m and RMS "
            f"{row['rms_residual_m']:.1f} m, against {raw_row['bias_m']:+.1f} m and "
            f"{raw_row['rms_residual_m']:.1f} m unsmoothed; it ranges "
            f"{row['residual_min_m']:+.1f} to {row['residual_max_m']:+.1f} m. What "
            f"survives this window is a place the {rs}–{re_} rate genuinely fails, at "
            "the scale the model resolves; what disappears between windows was "
            "transect-scale scatter. The correlation is deliberately not quoted: a "
            "symmetric smoother removes variance that is uncorrelated between the two "
            "sides, so r climbs with the window whether or not the smoothing is "
            "telling the truth (it is in tables/residual_by_scale.csv, named for that). "
            f"Seaward positive, in metres. {rf._marks_clause(s, e)}"
            f" The y axis is ±{half:g} m, fixed, the same on every metre figure."
            + offaxis_clause(off, half)))
        written += save(fig, sm["out"] / f"{prod.key}_smoothed_{s}_{e}_w{w:02d}")
        plt.close(fig)
    return written


def smooth_provenance(r, sm) -> None:
    s, e, years = r["start"], r["end"], r["years"]
    prod, (rs, re_) = r["prod"], r["rate_window"]
    st = sm["stats"]
    head = (f"| LOESS window | n domains | bias (m) | RMS residual (m) | residual range (m) "
            f"| sd {prod.tok} (m) | sd observed (m) | sign agreement | r |")
    lines = [head, "|" + "---|" * 9]
    for _, x in st.iterrows():
        w = int(x["window_domains"])
        lines.append(
            f"| {'raw (none)' if w == 0 else f'{w} domains / {x.window_km:g} km'} "
            f"| {int(x.n_domains)} | {x.bias_m:+.1f} | {x.rms_residual_m:.1f} "
            f"| {x.residual_min_m:+.1f} to {x.residual_max_m:+.1f} | {x.sd_rate_m:.1f} "
            f"| {x.sd_observed_m:.1f} | {x.pct_sign_agreement:.0f}% "
            f"| {x.r_inflated_by_smoothing:.2f} |")
    (sm["out"] / "PROVENANCE.md").write_text("\n".join([
        f"# 3-rates/coastsat/{prod.key}/{s}_{e}/smoothed - provenance",
        "",
        f"Written {dt.datetime.now():%Y-%m-%d %H:%M} by "
        "scripts/input_prep/5-scr/3-rates/coastsat/total_change/coastsat_total_change.py "
        f"(--product {prod.key}), beside the raw comparison one level up.",
        "",
        f"**{prod.noun}** = the {rs}-{re_} LRR x {years} yr"
        + (". The rate is fitted on the window it is evaluated over, so nothing is "
           "extrapolated." if (rs, re_) == (s, e) else
           f", evaluated over {s}-{e} -- a window the rate was NOT fitted on. That is "
           "what makes it a projection."),
        "",
        "## Why",
        "",
        "The rate the model is graded against is not the raw rate: it is the raw "
        f"domain mean over GIS 1-{SPLICE_DOMAINS} and a 10-domain alongshore LOESS of the "
        "transect values beyond (`cascade_pipeline.coastsat_loess`, imported here "
        "rather than re-implemented). The raw projected-vs-observed comparison "
        "therefore tests a quantity nobody feeds the model. This one tests the target "
        "as it is applied.",
        "",
        "LOESS commutes with the x years multiply, so smoothing the RATE and smoothing "
        "the DISTANCE are the same operation; nothing here turns on the "
        "order. What matters is that both sides get the same pass at the same window, "
        f"including the GIS 1-{SPLICE_DOMAINS} splice, so no residual is a smoothed "
        "quantity minus an unsmoothed one.",
        "",
        "## The sweep",
        "",
        *lines,
        "",
        "**Read the bias and the RMS residual, not r.** A symmetric smoother strips "
        "high-frequency variance that is uncorrelated between the two sides, so r "
        "rises with the window whether or not the smoothing is right; the column is "
        "named `r_inflated_by_smoothing` in `tables/residual_by_scale.csv` for that "
        "reason. The bias is close to smoothing-invariant and is the honest summary "
        "of whether the trend over- or under-predicts net change.",
        "",
        "A departure that survives the 10-domain window is a place the "
        f"{rs}-{re_} rate genuinely fails at the scale the model resolves. A departure "
        "that collapses between 3 and 10 domains was transect-scale scatter in the "
        "LRR estimate, the observed endpoint, or both.",
        "",
        "## The figures",
        "",
        f"`{prod.key}_smoothed_{s}_{e}_w<NN>.png` is one panel per window, the "
        "rate-derived distance against observed, both through that window's pass -- "
        "the per-place question, does the trend hold HERE.",
        "",
        f"`{prod.key}_smoothed_{s}_{e}_overlay.png` (Hannah, 2026-09-21) puts every "
        "LOESS width's curve on one panel and drops the observed side, which answers "
        "the other question: what the window does to the target. It is the rate figure "
        f"`3-rates/coastsat/lrr/{rs}_{re_}/smoothing_windows_{rs}_{re_}.png` in metres -- the "
        f"LOESS commutes with the x {years} yr multiply, so the curves have the same "
        "shape and only the units differ. It is drawn because metres is the unit the "
        "model and the dune line are read in, not because it is a different field.",
        "",
        "## Caveats",
        "",
        f"- The observed side is the thinner estimate: the LRR is fitted through a "
        f"median {r['t']['n_obs'].median():.0f} satellite positions per transect, while the "
        f"endpoint uses {r['t']['n_start'].median():.0f} positions in {s} and "
        f"{r['t']['n_end'].median():.0f} in {e}. Most of the noise the LOESS is "
        "removing is probably observed-side, which is why both sides are smoothed.",
        "- A 10-domain window is 5 km and the beach fills inside this record are 3-5 km "
        "wide (2014 at GIS 84-89, 2022 at GIS 6-15, 2022 at GIS 21-28). The fill "
        "signature in the residual is smeared at 10 domains and is clearest in the "
        "3-domain panel; do not read its disappearance at 10 as evidence it was noise.",
        f"- GIS 1-{SPLICE_DOMAINS} are unsmoothed at every window, so the largest positive "
        "residual in the raw product (GIS 1) is carried through unchanged by "
        "construction.",
        "",
    ]), encoding="utf-8")


def provenance(r) -> None:
    s, e, years, t, dom = r["start"], r["end"], r["years"], r["t"], r["dom"]
    prod, (rs, re_) = r["prod"], r["rate_window"]
    ok = dom.dropna(subset=["mean_observed_change_m"])
    diff = ok["mean_observed_minus_rate_m"]
    corr = np.corrcoef(ok["mean_rate_change_m"], ok["mean_observed_change_m"])[0, 1]
    no_obs = int(t["observed_change_m"].isna().sum())
    same = (rs, re_) == (s, e)
    (r["out"] / "PROVENANCE.md").write_text(NL.join([
        f"# 3-rates/coastsat/{prod.key}/{s}_{e} - provenance",
        "",
        f"Written {dt.datetime.now():%Y-%m-%d %H:%M} by "
        "scripts/input_prep/5-scr/3-rates/coastsat/total_change/coastsat_total_change.py "
        f"(--product {prod.key}).",
        "",
        "## Which product this is",
        "",
        (f"**Total shoreline change.** The rate is fitted on {rs}-{re_} and evaluated "
         f"over {s}-{e} -- the SAME window -- so nothing is extrapolated and this is "
         "not a projection. "
         + (f"`../../projected/{s}_{e}/` is the other reading of this window: the "
            f"{PROJECTED_RATE_WINDOW[0]}-{PROJECTED_RATE_WINDOW[1]} rate carried onto "
            "it instead of its own."
            if (s, e) in PROJECTED.windows else
            "There is no projected counterpart of this window: the rate IS the "
            f"{PROJECTED_RATE_WINDOW[0]}-{PROJECTED_RATE_WINDOW[1]} one, so a "
            "projection of it onto itself would be these same numbers.")
         if same else
         f"**Projected shoreline change.** The rate is fitted on {rs}-{re_} but "
         f"evaluated over {s}-{e}, a window it was NOT fitted on. That is what makes "
         "it a projection: the long-term trend asked what this window should have "
         f"looked like. The same window's OWN rate is in `../../total_change/{s}_{e}/`, "
         "and the difference between the two is how much the long-term rate misses "
         "this period by before the model is involved."),
        "",
        f"**{prod.noun}** = the transect's {rs}-{re_} LRR "
        f"(`../../lrr/{rs}_{re_}/transect_lrr_full.csv`) x {years} yr ({e} - {s}). Per "
        "domain the mean over its transects; it equals the LRR table's `mean_lrr` x "
        f"{years} to {r['check']:.1e} m/yr.",
        "",
        f"**Observed** = mean CoastSat position over calendar {e} minus the mean over "
        f"calendar {s}, per transect, from the same time series the LRR is fitted to. "
        "No rate anywhere in it. "
        f"Median positions per transect: {t['n_start'].median():.0f} in {s}, "
        f"{t['n_end'].median():.0f} in {e}. {no_obs} of {len(t)} transects have no "
        "position in one of the two years and no observed change.",
        "",
        f"Seaward positive. `observed_minus_{prod.tok}_m` > 0 means the shoreline "
        "ended more seaward than the trend implies.",
        "",
        "## Island summary",
        "",
        f"Domain mean {prod.tok} {dom['mean_rate_change_m'].mean():+.1f} m, observed "
        f"{ok['mean_observed_change_m'].mean():+.1f} m. Landward in "
        f"{int((dom['mean_rate_change_m'] < 0).sum())} of {N_DOMAINS} domains, observed "
        f"in {int((ok['mean_observed_change_m'] < 0).sum())} of {len(ok)}. Observed minus "
        f"{prod.tok} per domain: mean {diff.mean():+.1f} m, range {diff.min():+.1f} to "
        f"{diff.max():+.1f} m; r = {corr:.2f}.",
        "",
    ]), encoding="utf-8")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="The CoastSat LRR as a distance, beside the observed change. "
                    "TOTAL = the rate over the window it was fitted on; "
                    "PROJECTED = the 1996-2024 rate carried onto a window it was not.")
    ap.add_argument("--product", choices=["total", "projected", "both"], default="total",
                    help="total (default) -> 3-rates/coastsat/total_change/; "
                         "projected -> 3-rates/coastsat/projected/.")
    ap.add_argument("--windows", nargs="+", metavar="START_END")
    ap.add_argument("--smooth-windows", nargs="+", type=int, default=list(SMOOTH_WINDOWS),
                    metavar="N", help="LOESS window widths in domain units (1 = 500 m).")
    ap.add_argument("--no-smoothed", action="store_true",
                    help="Build the raw comparison only, skipping smoothed/.")
    a = ap.parse_args(argv)
    prods = [TOTAL, PROJECTED] if a.product == "both" else [PRODUCTS[
        "total_change" if a.product == "total" else "projected"]]
    apply_style()
    cache: dict = {}          # shared across products: the chainage is the same
    for prod in prods:
        asked = ([tuple(int(x) for x in w.split("_")) for w in a.windows]
                 if a.windows else list(prod.windows))
        # A window a product is not defined for is skipped loudly rather than
        # built: `projected/1996_2024` would be `total_change/1996_2024` under
        # another name, and two folders of identical numbers is the exact
        # confusion this rename was done to end.
        wins = [w for w in asked if w in prod.windows]
        for w in asked:
            if w not in prod.windows:
                print(f"skip  {prod.key} {w[0]}_{w[1]}: not a {prod.key} window "
                      f"(its rate window IS {w[0]}-{w[1]}, so that is total change)")
        print(f"== {prod.noun} -> {prod.root.relative_to(_REPO)}")
        # The overlays share one y bound across every window built for this
        # product, so they are drawn in a second pass once all series exist.
        built: list = []
        for s, e in wins:
            r = build(s, e, cache, prod)
            figure(r)
            provenance(r)
            d = r["dom"]
            rs, re_ = r["rate_window"]
            print(f"{s}_{e}  LRR {rs}-{re_} x{r['years']} yr  "
                  f"{prod.tok} {d['mean_rate_change_m'].mean():+6.1f} m  "
                  f"observed {d['mean_observed_change_m'].mean():+6.1f} m  "
                  f"obs-{prod.tok[:4]} {d['mean_observed_minus_rate_m'].mean():+5.1f} m  "
                  f"(check {r['check']:.1e})  -> {r['out'].relative_to(_REPO)}")
            if a.no_smoothed:
                continue
            sm = smooth(r, tuple(a.smooth_windows))
            smooth_figures(r, sm)
            smooth_provenance(r, sm)
            built.append((r, sm))
            for _, x in sm["stats"].iterrows():
                w = int(x["window_domains"])
                print(f"    LOESS {'raw   ' if w == 0 else f'{w:>2d} dom':<7}"
                      f"{'' if w == 0 else f'{x.window_km:>4.1f} km'}  "
                      f"bias {x.bias_m:+6.1f} m  RMS {x.rms_residual_m:5.1f} m  "
                      f"range {x.residual_min_m:+6.1f} to {x.residual_max_m:+6.1f} m  "
                      f"sign {x.pct_sign_agreement:3.0f}%  (r {x.r_inflated_by_smoothing:.2f})")
        if built:
            half, tick = overlay_bounds([sm for _, sm in built])
            for r, sm in built:
                for pth in smooth_overlay_figure(r, sm, half, tick):
                    print(f"    overlay  {WINDOW_ROLE.get((r['start'], r['end']), '?'):<24}"
                          f"+/-{half:g} m  -> {pth.relative_to(_REPO)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
