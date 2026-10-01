"""
The two smoothed halves-overlay sheets with the zeroBE CASCADE hindcast drawn over them.

    python scripts/analyze_output/compare_runs/hindcast_vs_observed/smoothed_lowess7_with_cascade.py [--window 7]

Imports the observed curves from 5-scr's smoothed_lowess7_vs_duneline.py and
adds the model line; writes output/comparisons/target_comparison/smoothed_lowess7_with_cascade/.
Details: README.md beside this script.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
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
sys.path.insert(0, str(Path(__file__).resolve().parent))

import rates_figures as rf  # noqa: E402
import smoothed_lowess7_vs_duneline as sl7  # noqa: E402
import total_change_vs_duneline as tcd  # noqa: E402
import rate_windows as rw  # noqa: E402
import target_comparison as tc  # noqa: E402
from rates_figures import cw, plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import MultipleLocator  # noqa: E402
from cascade_pipeline.coastsat_lowess import spliced_lowess_series  # noqa: E402
from cascade_pipeline.domains import DEFAULT_DOMAINS  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    COMPARISONS_ROOT, DOMAIN_AXIS_LABEL, INK, INK_MUTED, _title, apply_style,
    caption, compare_header, figsize, mark_offaxis, offaxis_clause, save,
    support_dir,
)

# --- CONFIG ------------------------------------------------------------------
OUT_ROOT = COMPARISONS_ROOT / "target_comparison" / "smoothed_lowess7_with_cascade"
# Dark green: distinct from the sign fill, the dune line and the raw dots, in greyscale too
C_MODEL = "#1b6b3a"
LW_MODEL = 1.5
HALVES = sl7.HALVES
# -----------------------------------------------------------------------------


# The zeroBE run's net change in metres, raw and smoothed like the observations
def model_series(window, years, lowess_window):
    spec = tc.UNSOLVED_RUNS[window]
    df, row = rw.load_model(window, spec, "coastsat", preset=tc.UNSOLVED_PRESET)
    if df is None:
        raise SystemExit(f"no run for {window}: {spec[0]}")
    if "zeroBE" not in str(row.get("run_name", "")):
        raise SystemExit(
            f"{row.get('run_name')} is not a zeroBE run; the caption claims no "
            "source/sink term was fitted anywhere, which would be false")
    idx = pd.RangeIndex(1, sl7.N + 1, name="domain_number")
    d = df.set_index("domain_number")["change_rate_m_yr"].reindex(idx) * years
    ids = np.asarray(idx, dtype=int)
    centres = (ids - 0.5) * DEFAULT_DOMAINS.domain_spacing_m
    sm, _ = spliced_lowess_series(ids, centres, d.to_numpy(float),
                                 lowess_window, skip=sl7.SPLICE)
    return d, sm.reindex(idx), row


# Interior mean of model minus target
def _bias(model, target):
    lo, hi = rw.INTERIOR
    m = model.loc[lo:hi].to_numpy(float)
    t = target.loc[lo:hi].to_numpy(float)
    ok = np.isfinite(m) & np.isfinite(t)
    return float((m[ok] - t[ok]).mean())


# One sheet: both halves, the observed pair from sl7 plus the model line
def figure(product, stem, what, data, models, lowess_window, out_dir):
    half, tick = tcd.Y_HALF_M, tcd.Y_TICK_M
    km = lowess_window * DEFAULT_DOMAINS.domain_spacing_m / 1000.0
    fig, axes = plt.subplots(2, 1, sharex=True, sharey=True,
                             constrained_layout=True,
                             figsize=figsize("double", height=5.6))
    off, rows = [], []
    for i, (ax, w) in enumerate(zip(axes, HALVES)):
        r = data[w]
        off += [(f"{w[0]}–{w[1]} {lab}", pts)
                for lab, pts in sl7._panel(ax, r, half, w, label=(i == 0))]
        sm_m = models[w]["smoothed"].to_numpy(float)
        x = np.asarray(models[w]["smoothed"].index, dtype=float)
        ax.plot(x, sm_m, color=C_MODEL, lw=LW_MODEL, zorder=13,
                solid_capstyle="round")
        off.append((f"{w[0]}–{w[1]} the smoothed CASCADE run",
                    mark_offaxis(ax, x, sm_m, half, color=C_MODEL)))
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
        Line2D([], [], color=INK, lw=tcd.LW),
        Line2D([], [], color=C_MODEL, lw=LW_MODEL),
        Line2D([], [], color=INK_MUTED, marker="o", ms=2.2, lw=0),
    ], [f"Shoreline target, LOWESS {km:g} km",
        f"Dune line target, LOWESS {km:g} km",
        "CASCADE (zeroBE, full management, groin off)",
        "Raw domain means (both targets)"], ncol=2)
    compare_header(fig, [
        f"All three curves LOWESS-smoothed at {lowess_window} domains ({km:g} km); "
        f"GIS 1–{sl7.SPLICE} kept raw",
        f"{what}   ·   dune line always the sub-period's own measured change"])

    per = "; ".join(
        f"{w[0]}–{w[1]} model − shoreline {_bias(models[w]['smoothed'], r['sm_shore']):+.1f} m, "
        f"model − dune line {_bias(models[w]['smoothed'], r['sm_dune']):+.1f} m"
        for w, r in rows)
    caption(fig, (
        "The two candidate targets and the CASCADE hindcast by GIS domain (1 at "
        "Cape Point, 90 at Pea Island), 1996–2010 above 2010–2024, all three "
        f"passed through the same alongshore LOWESS at {lowess_window} domains "
        f"({km:g} km). {what[0].upper() + what[1:]}; the dune line is always that "
        "half's own measured net change. **Dark green is CASCADE**: the zeroBE "
        "run of the matrix cell — full management, groin off — which carries NO "
        "source/sink term in ANY domain, the two ends included, so all 90 "
        "domains are the model's own response and neither target was fitted "
        "anywhere in it. That is what makes it readable against both. The "
        "coloured line and fill are the shoreline target, the black line the "
        "dune-line target, and the faint dots behind them the raw domain means. "
        f"GIS 1–{sl7.SPLICE} keep their raw values on every curve — the Oregon "
        "Inlet boundary treatment — so nothing there is the smoother's doing. "
        "The observed targets are smoothed at transect resolution and the model "
        "only exists per domain, so its LOWESS runs over 90 points rather than "
        "~900; at 3.5 km that barely changes the fitted curve, but it is not "
        f"literally the same operation. Interior means: {per}. Seaward positive, "
        f"in metres; the y axis is ±{half:g} m. "
        + rf._marks_clause(*HALVES[1]) + offaxis_clause(off, half)))
    written = save(fig, out_dir / f"{stem}_with_cascade_1996_2010_2024")
    plt.close(fig)
    return written, rows


# Run: model series, one sheet per shoreline reading, tables and PROVENANCE.md
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--window", type=int, default=sl7.WINDOW,
                    help="LOWESS width in domain units (1 = 500 m)")
    a = ap.parse_args(argv)
    apply_style()
    OUT_ROOT.mkdir(parents=True, exist_ok=True)

    # The model line for each half, and which runs it came from
    models, prov = {}, []
    for w in HALVES:
        raw, sm, row = model_series(w, float(w[1] - w[0]), a.window)
        models[w] = {"raw": raw, "smoothed": sm}
        prov.append(row)
    pd.DataFrame(prov).to_csv(OUT_ROOT / "runs_used.csv", index=False)

    # Both sheets, with every plotted value kept
    written, long, summary = [], [], []
    for product, stem, what in sl7.SHEETS:
        data = sl7.build(product, a.window)
        w_, rows = figure(product, stem, what, data, models, a.window, OUT_ROOT)
        written += w_
        for w, r in rows:
            long.append(pd.DataFrame({
                "product": product, "window": f"{w[0]}_{w[1]}",
                "lowess_domains": a.window,
                "domain_number": np.asarray(r["sm_shore"].index, dtype=int),
                "shoreline_target_smoothed_m": r["sm_shore"].to_numpy(float),
                "duneline_target_smoothed_m": r["sm_dune"].to_numpy(float),
                "cascade_raw_m": models[w]["raw"].to_numpy(float),
                "cascade_smoothed_m": models[w]["smoothed"].to_numpy(float),
            }))
            summary.append({
                "product": product, "window": f"{w[0]}_{w[1]}",
                "bias_vs_shoreline_m": round(_bias(models[w]["smoothed"], r["sm_shore"]), 1),
                "bias_vs_duneline_m": round(_bias(models[w]["smoothed"], r["sm_dune"]), 1),
            })
    pd.concat(long, ignore_index=True).round(3).to_csv(
        OUT_ROOT / "domain_values.csv", index=False)
    summ = pd.DataFrame(summary)
    summ.to_csv(support_dir(OUT_ROOT) / "island_summary.csv", index=False)

    # PROVENANCE.md with the interior means
    tbl = ["| shoreline reading | window | model − shoreline (m) | model − dune line (m) |",
           "|" + "---|" * 4]
    for _, x in summ.iterrows():
        tbl.append(f"| {x['product']} | {x['window'].replace('_', '–')} "
                   f"| {x.bias_vs_shoreline_m:+.1f} | {x.bias_vs_duneline_m:+.1f} |")
    (OUT_ROOT / "PROVENANCE.md").write_text("\n".join([
        "# smoothed_lowess7_with_cascade — provenance", "",
        f"Written {dt.datetime.now():%Y-%m-%d %H:%M} by "
        "`scripts/analyze_output/compare_runs/hindcast_vs_observed/smoothed_lowess7_with_cascade.py`.",
        "",
        "The observations-only pair in "
        "`data/hatteras_init/5-scr/4-comparisons/shoreline_vs_duneline/"
        "smoothed_lowess7/` with the CASCADE hindcast over it in dark green. "
        "Every detail of the two observed curves is imported from that script, "
        "so the sheets differ in exactly one thing.", "",
        "## The run", "",
        "zeroBE, full management, groin off, one per period — "
        "`HAT_1996_2010_zeroBE_offsetmetres_road_bdm_nogroin` and "
        "`HAT_2010_2024_zeroBE_offsetmetres_road_bdm_nourish_nogroin` (see runs_used.csv). "
        "NO source/sink term in any domain, the two ends included, so all 90 "
        "are the model's own response and NEITHER target was fitted anywhere "
        "in it. Chosen over the headline edgeBE matrix run for exactly that "
        "reason: edgeBE solves GIS 1 and 90 against the CoastSat target, which "
        "would make the model partly fitted to one of the two things it is "
        "being compared against.", "",
        "## Interior means, GIS 2–89", "",
        *tbl, "",
        "## Reading it", "",
        f"All three curves are smoothed at {a.window} domains "
        f"({a.window * DEFAULT_DOMAINS.domain_spacing_m / 1000.0:g} km) so no "
        "pair is a smoothed quantity against an unsmoothed one. The observed "
        "sides are smoothed at TRANSECT resolution and the model only exists "
        "per domain, so its LOWESS runs over 90 points rather than ~900 — at "
        "this width that barely changes the fitted curve, but the two are not "
        "literally the same operation.", "",
        f"GIS 1–{sl7.SPLICE} are unsmoothed on every curve, so no difference "
        "there is the smoother's.", "",
        "Correlations are deliberately not quoted: a symmetric smoother "
        "inflates r on both sides regardless of whether the curves genuinely "
        "agree better at this width, which "
        "`../smoothing_scale/PROVENANCE.md` established over 96 rows. Read "
        "the bias.", "",
    ]), encoding="utf-8")
    for p in written:
        print("wrote    ", p.relative_to(_REPO))
    print()
    print(summ.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
