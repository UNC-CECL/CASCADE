"""
The CoastSat LRR turned into a distance, beside the distance the shoreline actually moved.

    python scripts/input_prep/5-scr/3-rates/coastsat/total_change/coastsat_total_change.py
    python scripts/input_prep/5-scr/3-rates/coastsat/total_change/coastsat_total_change.py --product projected
    python scripts/input_prep/5-scr/3-rates/coastsat/total_change/coastsat_total_change.py --windows 1996_2024 2010_2024

TOTAL uses the window's own rate; PROJECTED the 1996-2024 rate. Writes
figures, tables and provenance per window, raw and smoothed. Details: scripts/input_prep/5-scr/3-rates/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
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
# The model target's own smoother, so this figure is the graded treatment
from cascade_pipeline.coastsat_lowess import spliced_lowess_series  # noqa: E402
from cascade_pipeline.domains import DEFAULT_DOMAINS  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
N_DOMAINS = 90
# The canonical 1996 -> 2010 -> 2024 chain
OBS_LW = 1.1
# The observed CoastSat line, PURPLE since 2026-09-22 (Hannah)
C_OBSERVED = C["ACCENT"]
# ONE fixed metre axis on every figure here (Hannah, 2026-09-22)
Y_HALF_M = 100.0
Y_TICK_M = 20.0
NL = chr(10)          # the provenance writers join on it
# -----------------------------------------------------------------------------


# One of the two named products
class Product:

    def __init__(self, key, noun, tok, root, windows, rate_window, method):
        self.key, self.noun, self.tok = key, noun, tok
        self.root, self.windows = root, windows
        self._rate_window, self._method = rate_window, method

    # The window the LRR is FITTED on, which is what names the product
    def rate_window(self, s, e):
        return self._rate_window(s, e)

    # The parenthetical in a figure title: where the rate came from and how the distance was made
    def method(self, s, e):
        return self._method(s, e)

    # Internal `rate` column names -> the product's written ones
    def cols(self, df):
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
# The 1996-2024 rate carried onto a window it was not fitted on
PROJECTED = Product(
    "projected", "Projected shoreline change", "projected", COASTSAT_PROJECTED_ROOT,
    [(1996, 2010), (2010, 2024)],
    rate_window=lambda s, e: PROJECTED_RATE_WINDOW,
    method=lambda s, e: (f"CoastSat LRR {PROJECTED_RATE_WINDOW[0]}–"
                         f"{PROJECTED_RATE_WINDOW[1]} × {e - s} yr"),
)
PRODUCTS = {p.key: p for p in (TOTAL, PROJECTED)}

# LOWESS window widths in domain units (1 domain = 500 m)
SMOOTH_WINDOWS = (3, 5, 7, 10)
# GIS 1..SPLICE_DOMAINS keep their raw domain means instead of the LOWESS
SPLICE_DOMAINS = 10


# (mean position, n) over one calendar year, or (nan, 0)
def year_mean(df, year):
    sel = df.loc[df["date"].dt.year == year, "chainage"]
    return (float(sel.mean()), int(sel.size)) if sel.size else (np.nan, 0)


# The product's rate as a distance over start-end, beside the observed change over the same years
def build(start: int, end: int, cache: dict, prod: Product = TOTAL) -> dict:
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
    # Internal names are neutral (`rate`)
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


# One window's rate-as-distance against the observed change
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
    # Quantity, window, method (Hannah, 2026-09-21)
    ax.set_title(f"{prod.noun}, {s}–{e} ({prod.method(s, e)})",
                 pad=20 if cw.fills_in(s, e) else 6)
    # Both series here are CoastSat; the observed line is purple
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


# One alongshore LOWESS pass at transect resolution, averaged to domains, GIS 1..SPLICE_DOMAINS raw
def _smooth_series(dom_ids, along_m, values, window):
    return spliced_lowess_series(dom_ids, along_m, values, window,
                                skip=SPLICE_DOMAINS)


# The rate-vs-observed comparison under the target's LOWESS at each window; window 0 is raw
def smooth(r, windows=SMOOTH_WINDOWS) -> dict:
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
            lowess_frac=frac, n_domains=int(ok.sum()),
            bias_m=float(d.mean()), rms_residual_m=float(np.sqrt((d ** 2).mean())),
            residual_min_m=float(d.min()), residual_max_m=float(d.max()),
            sd_rate_m=float(a.std(ddof=1)), sd_observed_m=float(b.std(ddof=1)),
            pct_sign_agreement=float(100.0 * (np.sign(a) == np.sign(b)).mean()),
            # Named for what it is: smoothing inflates r, so it is not a score
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


# The fixed metre axis, as everywhere else here (Hannah, 2026-09-22)
def _smooth_bounds(sm):
    return Y_HALF_M, Y_TICK_M


# One y half-range and tick for EVERY window's overlay, from the rate-derived series alone
def overlay_bounds(sms):
    return Y_HALF_M, Y_TICK_M


# Every LOWESS width's distance on ONE panel, no observed side (Hannah, 2026-09-21)
def smooth_overlay_figure(r, sm, half=None, tick=None) -> list:
    s, e, years = r["start"], r["end"], r["years"]
    prod, (rs, re_) = r["prod"], r["rate_window"]
    if half is None:
        half, tick = overlay_bounds([sm])
    windows = (0,) + tuple(sm["windows"])
    km_of = DEFAULT_DOMAINS.domain_spacing_m / 1000.0
    role = WINDOW_ROLE.get((s, e))

    fig, ax = plt.subplots(figsize=figsize("double", aspect=0.40),
                           constrained_layout=True)
    # Mean all-NaN draws the frame, grid, village bands and structures with no sign fill
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
    # Quantity, window, method, like every other figure in the tree since 2026-09-21
    ax.set_title(f"{prod.noun}, {s}–{e} ({prod.method(s, e)}, every LOWESS width)",
                 pad=20 if fills else 6)

    h = [Line2D([], [], color=SMOOTH_RAMP[0], lw=0.8, marker="o", ms=2.2)]
    h += [Line2D([], [], color=SMOOTH_RAMP[i % len(SMOOTH_RAMP)], lw=1.5)
          for i, w in enumerate(windows) if w]
    labels = ["Unsmoothed domain means"] + [
        f"LOWESS {w * km_of:g} km ({w} domains)" for w in windows if w]
    # At most three legend columns
    fig.legend(h, labels, loc="outside lower center", ncol=min(len(h), 3), frameon=False)

    # How far apart the windows are, in the unit of the axis.
    spread = {w: float((sm["series"][w][0] - sm["series"][0][0]).abs().max())
              for w in windows if w}
    spread_txt = "; ".join(f"{w * km_of:g} km up to {v:.0f} m" for w, v in spread.items())
    # The rate figure exists only where coastsat_lrr_smoothing_windows.py has been run
    rate_png = (COASTSAT_LRR_ROOT / f"{rs}_{re_}" / f"smoothing_windows_{rs}_{re_}.png")
    rate_ref = (
        f"This is the rate figure `3-rates/coastsat/lrr/{rs}_{re_}/{rate_png.name}` in "
        f"metres — a LOWESS commutes with the × {years} yr multiply, so the curves have "
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
        f"LOWESS the model target is built through, at {', '.join(f'{w * km_of:g} km' for w in windows[1:-1])} "
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


# One panel per window, all on the same y axis so they read side by side
def smooth_figures(r, sm) -> list:
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
        labels = [f"{prod.noun}, LOWESS {km:g} km ({prod.method(s, e)})",
                  f"{prod.noun}, individual transects (raw)",
                  f"CoastSat observed change, LOWESS {km:g} km"]
        fig.legend(h, labels, loc="outside lower center", ncol=2, frameon=False,
                   handler_map={tuple: HandlerTuple(ndivide=None, pad=0.3)})
        # Quantity, window, method -- plus the smoothing width
        ax.set_title(f"{prod.noun}, {s}–{e} ({prod.method(s, e)}, LOWESS {km:g} km)",
                     pad=20 if cw.fills_in(s, e) else 6)
        # Both series here are CoastSat; the observed line is purple
        compare_header(fig, [
            f"{s}–{e}   ·   BOTH series are CoastSat — no dune line on this figure",
            f"{prod.noun.lower()}: {prod.method(s, e)}   ·   observed: CoastSat mean position, all of {e} minus all of {s}"])
        caption(fig, (
            f"{prod.noun} and observed change, {s}–{e}, both "
            f"passed through the same alongshore LOWESS of {w} domains ({km:g} km). "
            "The quantity the model is graded against is not the raw rate but this "
            "one — raw domain means over GIS 1–10 and a LOWESS of the transect values "
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


# The smoothed comparison's provenance
def smooth_provenance(r, sm) -> None:
    s, e, years = r["start"], r["end"], r["years"]
    prod, (rs, re_) = r["prod"], r["rate_window"]
    st = sm["stats"]
    head = (f"| LOWESS window | n domains | bias (m) | RMS residual (m) | residual range (m) "
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
        f"domain mean over GIS 1-{SPLICE_DOMAINS} and a 10-domain alongshore LOWESS of the "
        "transect values beyond (`cascade_pipeline.coastsat_lowess`, imported here "
        "rather than re-implemented). The raw projected-vs-observed comparison "
        "therefore tests a quantity nobody feeds the model. This one tests the target "
        "as it is applied.",
        "",
        "LOWESS commutes with the x years multiply, so smoothing the RATE and smoothing "
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
        "LOWESS width's curve on one panel and drops the observed side, which answers "
        "the other question: what the window does to the target. It is the rate figure "
        f"`3-rates/coastsat/lrr/{rs}_{re_}/smoothing_windows_{rs}_{re_}.png` in metres -- the "
        f"LOWESS commutes with the x {years} yr multiply, so the curves have the same "
        "shape and only the units differ. It is drawn because metres is the unit the "
        "model and the dune line are read in, not because it is a different field.",
        "",
        "## Caveats",
        "",
        f"- The observed side is the thinner estimate: the LRR is fitted through a "
        f"median {r['t']['n_obs'].median():.0f} satellite positions per transect, while the "
        f"endpoint uses {r['t']['n_start'].median():.0f} positions in {s} and "
        f"{r['t']['n_end'].median():.0f} in {e}. Most of the noise the LOWESS is "
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


# What the figure's numbers are and where they came from
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


# Run: every window and product
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
                    metavar="N", help="LOWESS window widths in domain units (1 = 500 m).")
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
        # A window a product is not defined for is skipped loudly rather than built
        wins = [w for w in asked if w in prod.windows]
        for w in asked:
            if w not in prod.windows:
                print(f"skip  {prod.key} {w[0]}_{w[1]}: not a {prod.key} window "
                      f"(its rate window IS {w[0]}-{w[1]}, so that is total change)")
        print(f"== {prod.noun} -> {prod.root.relative_to(_REPO)}")
        # The overlays share one y bound across every window built for this product
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
                print(f"    LOWESS {'raw   ' if w == 0 else f'{w:>2d} dom':<7}"
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
