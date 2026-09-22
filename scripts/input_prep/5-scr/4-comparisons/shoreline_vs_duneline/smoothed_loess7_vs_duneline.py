"""
smoothed_loess7_vs_duneline.py
==============================================================================
The two halves-overlay sheets again, with BOTH curves passed through the
model target's alongshore LOESS at a 7-domain (3.5 km) window. Built
2026-09-22 (Hannah, by interview) as a place to see what the smoother does to
the shoreline-vs-dune-line comparison, in the two readings of the shoreline
side side by side.

WHAT IS DRAWN, one sheet per shoreline reading, 1996-2010 above 2010-2024:

    projected     shoreline = the 1996-2024 LRR x 14 yr, the SAME in both
                  panels (the long-term trend carried onto each half)
    total_change  shoreline = each half's OWN LRR x its own 14 yr

    The dune side is the same either way and always follows the sub-period:
    the measured net change between the two digitized lines bounding that
    half. Only the shoreline reading differs between the two sheets.

THE SMOOTHING (Hannah's choices, 2026-09-22)

    window        7 domains = 3.5 km.
    both sides    BOTH curves get the same pass at the same window, at
                  TRANSECT resolution, then average to domains. Smoothing one
                  side and not the other would make the gap between them an
                  artefact of the treatment rather than a beach-width change
                  -- the same trap the 3-rates smoothed panels avoid.
                  CoastSat has ~10 transects per domain and the dune line
                  exactly 5; `rates_figures._along` gives both an even spread
                  inside their domain, so the two are handled identically.
    GIS 1-10      kept at their RAW domain means, the Oregon Inlet boundary
                  treatment (`coastsat_loess.LoessConfig.skip_southern_domains`).
                  Hannah chose to keep it so the figure shows the target the
                  way the model actually sees it. The cost, stated here so it
                  is not read as a result: those ten domains are IDENTICAL to
                  the unsmoothed sheet by construction, and any difference
                  there is not the smoother.
    raw kept      the raw domain means stay on the figure as faint dots
                  behind each curve, so what the smoother removed is visible
                  without leaving the sheet.

WHAT THIS IS FOR.  It is a test folder, not a product: the question is how
much of the shoreline-vs-dune-line disagreement survives smoothing at the
scale the model resolves. The unsmoothed sheets it is paired with are
`coastsat_{projected,total_change}_vs_duneline_endpoint/all_windows_stacked/
*_halves_overlay.png`.

OUTPUT   data/hatteras_init/5-scr/4-comparisons/shoreline_vs_duneline/smoothed_loess7/
    loess7_projected_vs_duneline_1996_2010_2024_halves_overlay.png
    loess7_total_change_vs_duneline_1996_2010_2024_halves_overlay.png
    domain_smoothed.csv     per domain per window per product: both sides raw
                            and smoothed, and the beach-width gap of each
    PROVENANCE.md           what changed, per window and per product
    README.md, supporting/  (PDFs, CAPTIONS.md)

USAGE
    python scripts/input_prep/5-scr/4-comparisons/shoreline_vs_duneline/smoothed_loess7_vs_duneline.py
    python ... --window 7          # LOESS width in domain units
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

import rates_figures as rf  # noqa: E402
import total_change_vs_duneline as tcd  # noqa: E402
from rates_figures import cw, plt  # noqa: E402
from coastsat_vs_duneline import beach_width_handles, shade_beach_width  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import MultipleLocator  # noqa: E402
from cascade_pipeline.coastsat_loess import spliced_loess_series  # noqa: E402
from cascade_pipeline.domains import DEFAULT_DOMAINS  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    DOMAIN_AXIS_LABEL, INK, INK_MUTED, _title, apply_style, caption,
    compare_header, figsize, mark_offaxis, offaxis_clause, save, support_dir,
)
from site_layer.hat_observed_rates import (  # noqa: E402
    DUNELINE_ENDPOINT_ROOT, ENDPOINT_TRANSECT_FILE, SHORELINE_VS_DUNELINE,
)

OUT_ROOT = SHORELINE_VS_DUNELINE / "smoothed_loess7"
WINDOW = 7                      # domain units; 7 x 500 m = 3.5 km
SPLICE = 10                     # GIS 1..SPLICE keep their raw domain means
HALVES = [(1996, 2010), (2010, 2024)]
N = cw.N_DOMAINS
RAW_DOT_S = 3.0                 # the faint raw domain means behind each curve
RAW_ALPHA = 0.45
# Which shoreline readings to draw, and how each is said in one line.
SHEETS = [
    ("projected", "loess7_projected_vs_duneline",
     "shoreline is the SAME CoastSat LRR 1996–2024 × 14 yr in both panels"),
    ("total", "loess7_total_change_vs_duneline",
     "shoreline is each panel's OWN CoastSat LRR × its own 14 yr"),
]


def smooth_side(frame, value_col, window):
    """One alongshore LOESS pass at TRANSECT resolution, averaged to domains.

    The same two steps the scoring target is built through, so the curve here
    is the quantity the model is graded against rather than a different
    smoother that happens to look similar.
    """
    t, x_dom = rf._along(frame)
    along_m = x_dom * DEFAULT_DOMAINS.domain_spacing_m
    series, _ = spliced_loess_series(
        t["domain_number"].to_numpy(int), along_m,
        t[value_col].to_numpy(float), window, skip=SPLICE)
    return series.reindex(pd.RangeIndex(1, N + 1, name="domain_number"))


def build(product, window):
    """Per half: the raw domain means and the smoothed series, both sides."""
    tcd.PROD = tcd.PRODUCTS[product]
    out = {}
    for w in HALVES:
        dom, cs, meta, years, dune_years = tcd.load(w)
        dune_t = pd.read_csv(DUNELINE_ENDPOINT_ROOT / f"{w[0]}_{w[1]}"
                             / ENDPOINT_TRANSECT_FILE)
        dune_t = dune_t[dune_t["domain_number"].between(1, N)].copy()
        idx = pd.RangeIndex(1, N + 1, name="domain_number")
        d = dom.set_index("domain_number")
        out[w] = {
            "meta": meta, "years": years, "dune_years": dune_years,
            "raw_shore": d["shoreline_change_m"].reindex(idx),
            "raw_dune": d["dune_change_m"].reindex(idx),
            "sm_shore": smooth_side(cs, "shoreline_change_m", window),
            "sm_dune": smooth_side(dune_t, "change_m", window),
        }
    return out


def _panel(ax, r, half, window_tuple, label):
    """The smoothed pair with the raw domain means faint behind them."""
    sm_s = r["sm_shore"].to_numpy(float)
    sm_d = r["sm_dune"].to_numpy(float)
    x = np.asarray(r["sm_shore"].index, dtype=float)

    frame = pd.DataFrame({"domain_number": x.astype(int),
                          "mean_lrr": sm_s, "std_lrr": 0.0})
    cw.draw_panel(ax, frame, half, label=label, std=False)
    # raw first, so the smoothed curves sit over them
    raw_s = r["raw_shore"].to_numpy(float)
    ax.scatter(x, np.where(np.abs(raw_s) <= half, raw_s, np.nan), s=RAW_DOT_S,
               c=np.where(raw_s < 0, cw.C_ERODE, cw.C_ACCRETE),
               alpha=RAW_ALPHA, linewidths=0, zorder=3.5)
    raw_d = r["raw_dune"].to_numpy(float)
    ax.scatter(x, np.where(np.abs(raw_d) <= half, raw_d, np.nan), s=RAW_DOT_S,
               c=INK_MUTED, alpha=RAW_ALPHA, linewidths=0, zorder=3.6)
    ax.plot(x, sm_d, color=INK, lw=tcd.LW, zorder=12, solid_capstyle="round")
    cw.draw_shoals(ax, label=label)
    if cw.fills_in(*window_tuple):
        cw.draw_fills(ax, cw.fills_in(*window_tuple), half)
    return [("the smoothed shoreline", mark_offaxis(ax, x, sm_s, half, color=INK)),
            ("the smoothed dune line", mark_offaxis(ax, x, sm_d, half, color=INK))]


def figure(product, stem, what, data, window, out_dir):
    half, tick = tcd.Y_HALF_M, tcd.Y_TICK_M
    km = window * DEFAULT_DOMAINS.domain_spacing_m / 1000.0
    fig, axes = plt.subplots(2, 1, sharex=True, sharey=True,
                             constrained_layout=True,
                             figsize=figsize("double", height=5.6))
    off, rows = [], []
    for i, (ax, w) in enumerate(zip(axes, HALVES)):
        r = data[w]
        off += [(f"{w[0]}–{w[1]} {lab}", pts)
                for lab, pts in _panel(ax, r, half, w, label=(i == 0))]
        ax.yaxis.set_major_locator(MultipleLocator(tick))
        m = r["meta"]
        _title(ax, i, f"{w[0]}–{w[1]}   ·   dune line {m['start_date']} → "
                      f"{m['end_date']}"
                      + (" (assumed)" if bool(m["end_date_assumed"]) else "")
                      + f",  {r['dune_years']:.1f} yr")
        tcd._pad_title(ax, w)
        rows.append((w, r))
    axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)
    fig.supylabel(tcd.Y_LABEL, fontsize=9)

    tcd._legend(fig, [
        (Line2D([], [], color=cw.C_ACCRETE, lw=1.0),
         Line2D([], [], color=cw.C_ERODE, lw=1.0)),
        (Line2D([], [], color=cw.C_ACCRETE, marker="o", ms=2.2, lw=0),
         Line2D([], [], color=cw.C_ERODE, marker="o", ms=2.2, lw=0)),
        Line2D([], [], color=INK, lw=tcd.LW),
        Line2D([], [], color=INK_MUTED, marker="o", ms=2.2, lw=0),
    ], [f"Shoreline, LOESS {km:g} km", "Shoreline, raw domain means",
        f"Dune line, LOESS {km:g} km", "Dune line, raw domain means"], ncol=4)
    compare_header(fig, [
        f"BOTH curves LOESS-smoothed at {window} domains ({km:g} km); "
        f"GIS 1–{SPLICE} kept raw",
        f"{what}   ·   dune line always the sub-period's own measured change"])

    per = "; ".join(
        f"{w[0]}–{w[1]} r {_r(r['raw_shore'], r['raw_dune']):.2f} → "
        f"{_r(r['sm_shore'], r['sm_dune']):.2f}, mean beach width "
        f"{(r['raw_shore'] - r['raw_dune']).mean():+.1f} → "
        f"{(r['sm_shore'] - r['sm_dune']).mean():+.1f} m"
        for w, r in rows)
    caption(fig, (
        f"Shoreline against dune line by GIS domain (1 at Cape Point, 90 at Pea "
        f"Island), 1996–2010 above 2010–2024, with **both** curves passed through "
        f"the model target's alongshore LOESS at {window} domains ({km:g} km). "
        f"{what[0].upper() + what[1:]}; the dune line is always that half's own "
        "measured net change between its two digitized lines. The coloured line "
        "and fill are the smoothed shoreline (blue seaward, red landward), the "
        "black line the smoothed dune line, and the faint dots behind each are "
        "the RAW domain means, so what the smoother removed stays visible. Both "
        "sides get the same pass at the same width, at transect resolution, so "
        "the gap between them is not an artefact of treating them differently. "
        f"GIS 1–{SPLICE} keep their raw domain means — the boundary treatment at "
        "Oregon Inlet that the scoring target uses — so those ten domains are "
        "identical to the unsmoothed sheet by construction and nothing there is "
        f"the smoother's doing. Raw → smoothed, per window: {per}. Seaward "
        f"positive, in metres; the y axis is ±{half:g} m. "
        + rf._marks_clause(*HALVES[1]) + offaxis_clause(off, half)))
    written = save(fig, out_dir / f"{stem}_1996_2010_2024_halves_overlay")
    plt.close(fig)
    return written, rows


def _r(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    return float(np.corrcoef(a[ok], b[ok])[0, 1])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--window", type=int, default=WINDOW,
                    help="LOESS width in domain units (1 = 500 m)")
    a = ap.parse_args(argv)
    apply_style()
    OUT_ROOT.mkdir(parents=True, exist_ok=True)

    written, long, summary = [], [], []
    for product, stem, what in SHEETS:
        data = build(product, a.window)
        w_, rows = figure(product, stem, what, data, a.window, OUT_ROOT)
        written += w_
        for w, r in rows:
            long.append(pd.DataFrame({
                "product": product, "window": f"{w[0]}_{w[1]}",
                "loess_domains": a.window,
                "domain_number": np.asarray(r["sm_shore"].index, dtype=int),
                "shoreline_raw_m": r["raw_shore"].to_numpy(float),
                "shoreline_smoothed_m": r["sm_shore"].to_numpy(float),
                "duneline_raw_m": r["raw_dune"].to_numpy(float),
                "duneline_smoothed_m": r["sm_dune"].to_numpy(float),
                "beach_width_raw_m": (r["raw_shore"] - r["raw_dune"]).to_numpy(float),
                "beach_width_smoothed_m": (r["sm_shore"] - r["sm_dune"]).to_numpy(float),
            }))
            summary.append({
                "product": product, "window": f"{w[0]}_{w[1]}",
                "r_raw": round(_r(r["raw_shore"], r["raw_dune"]), 3),
                "r_smoothed": round(_r(r["sm_shore"], r["sm_dune"]), 3),
                "beach_width_raw_m": round(float((r["raw_shore"] - r["raw_dune"]).mean()), 1),
                "beach_width_smoothed_m": round(float((r["sm_shore"] - r["sm_dune"]).mean()), 1),
            })
    pd.concat(long, ignore_index=True).round(3).to_csv(
        OUT_ROOT / "domain_smoothed.csv", index=False)
    summ = pd.DataFrame(summary)
    summ.to_csv(support_dir(OUT_ROOT) / "island_summary.csv", index=False)

    tbl = ["| shoreline reading | window | r raw | r smoothed | beach width raw (m) "
           "| beach width smoothed (m) |", "|" + "---|" * 6]
    for _, x in summ.iterrows():
        tbl.append(f"| {x['product']} | {x['window'].replace('_', '–')} | {x.r_raw:.2f} "
                   f"| {x.r_smoothed:.2f} | {x.beach_width_raw_m:+.1f} "
                   f"| {x.beach_width_smoothed_m:+.1f} |")
    (OUT_ROOT / "PROVENANCE.md").write_text("\n".join([
        "# smoothed_loess7 — provenance", "",
        f"Written {dt.datetime.now():%Y-%m-%d %H:%M} by "
        "`scripts/input_prep/5-scr/4-comparisons/shoreline_vs_duneline/smoothed_loess7_vs_duneline.py`.", "",
        f"Both curves LOESS-smoothed at **{a.window} domains "
        f"({a.window * DEFAULT_DOMAINS.domain_spacing_m / 1000.0:g} km)**, at "
        f"transect resolution, with GIS 1–{SPLICE} kept at their raw domain "
        "means (the scoring target's Oregon Inlet treatment, Hannah's choice "
        "2026-09-22). The dune line always follows the sub-period.", "",
        *tbl, "",
        "**r rises with smoothing on both sides.** That is what a symmetric "
        "smoother does — it strips high-frequency variance that is uncorrelated "
        "between the two series — so it is NOT evidence that the two features "
        "agree better at 3.5 km. Read the beach width, which is close to "
        "smoothing-invariant, and compare the SHAPE against the unsmoothed "
        "sheets in "
        "`coastsat_{projected,total_change}_vs_duneline_endpoint/"
        "all_windows_stacked/`.", "",
        f"GIS 1–{SPLICE} are unsmoothed by construction, so no difference there "
        "is the smoother's.", "",
    ]), encoding="utf-8")
    for p in written:
        print("wrote    ", p.relative_to(_REPO))
    print()
    print(summ.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
