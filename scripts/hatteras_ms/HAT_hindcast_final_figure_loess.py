#!/usr/bin/env python3
"""The calibrated hindcast against the LOESS reference curve — presentation figure.

A companion to HAT_hindcast_final_figure.py, differing in three deliberate ways.

    SHARED Y AXIS
        Both periods on identical limits, so the panels can be read against
        each other rather than only against their own observations. That
        comparison is the point: period 2 sits almost entirely above period 1,
        which is the post-Isabel recovery and the nourishment era showing up as
        a whole-island shift in the rate, not as a local feature.

    THE LOESS CURVE ONLY, NOT THE SPLICED TARGET
        The calibration target is not one curve: GIS 1-10 are raw per-domain
        means and D11 north is the 10-domain LOESS. That splice is right for
        calibrating -- the raw means keep the short-wavelength signal the
        source/sink field has to answer for -- but it makes an awkward figure,
        because the eye reads a change of estimator as a change of coast.
        Here the LOESS is drawn throughout.

        D1-D10 IS DRAWN DASHED. The project excludes the LOESS there by
        convention (`skip_southern_domains = 10`), because the smoother is
        poorly constrained at the end of its range and Cape Point's
        attachment-detachment cycle is exactly the short-wavelength signal a
        5 km smoother destroys. Dashing it shows the data without implying it
        carries the same weight.

    SCORED ON D11-D89
        Where the LOESS curve and the calibration target are the SAME numbers,
        so the statistic on the panel is the statistic for the line drawn. The
        companion figure scores D2-D89 against the spliced target and reports
        0.526 / 0.558 m/yr; both are correct for what they measure, and mixing
        them would be the error.

ALSO ADDED FOR PRESENTATION
    An observational spread band (+/- 1 SD of transect rates within each
    domain, from `unc_m_yr`'s parent transect table). A domain is a 500 m
    average over ~10 transects, and showing that spread makes clear how much
    real alongshore variability the domain mean hides -- and therefore which
    model-observation gaps are meaningful and which sit inside the noise.

    Place names along the top, from PHYSICAL_ZONES, so an audience can locate
    features without a separate map.

Usage:
    python HAT_hindcast_final_figure_loess.py [--preset calibBE]

Writes output/comparisons/hindcast_calibrated_loess_reference.png

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import argparse
import importlib.util
import pathlib
import sys

import numpy as np
import pandas as pd

_HERE = pathlib.Path(__file__).resolve()
PROJECT_BASE_DIR = next(p for p in _HERE.parents if (p / "pyproject.toml").exists())
RAW_RUNS = PROJECT_BASE_DIR / "output" / "raw_runs"
OUT_DIR = PROJECT_BASE_DIR / "output" / "comparisons"

sys.path.insert(0, str(PROJECT_BASE_DIR / "scripts"))

from cascade_pipeline.run_registry import find_run_dir          # noqa: E402

LOESS_PATH = (PROJECT_BASE_DIR / "scripts" / "input_prep" / "7-source-sink"
              / "loess_smooth" / "HAT_be_zone_LOESS_analysis.py")

PERIODS = {
    "1984_2004": dict(panel="a", label="1984–2004", start=1984,
                      scenario="road_bdm", colour="#1565C0"),
    "2004_2024": dict(panel="b", label="2004–2024", start=2004,
                      scenario="road_bdm_nourish", colour="#B71C1C"),
}
LOCKED = (1, 90)
SKIP_SOUTH = 10               # matches LOESS_CONFIG.skip_southern_domains
SCORE_DOMAINS = range(SKIP_SOUTH + 1, 90)
RESERVED_COLOUR = "#FF8C00"

# Non-overlapping display spans. PHYSICAL_ZONES overlaps at D9-D10 (Cape Point
# and Buxton-Avon both claim them) and assign_physical_zone resolves that by
# first match; a label strip has to pick one, so it picks the same one.
PLACE_LABELS = [
    (1, 10, "Cape Point"), (11, 20, "Buxton"), (21, 31, "Avon"),
    (32, 59, "Mid-island"), (60, 74, "Wimble Shoals"),
    (75, 83, "Rodanthe"), (84, 90, "Pea Island"),
]


def analysis_module():
    spec = importlib.util.spec_from_file_location("_loess", LOESS_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def loess_and_spread(module, start, csv_path):
    """(loess_by_domain, sd_by_domain) for one period.

    The LOESS is the widest window's smoothed curve across ALL domains --
    build_target_table would splice raw means over D1-D10, which is what this
    figure is deliberately not doing.
    """
    from cascade_pipeline.coastsat_loess import (CoastSatDataset,
                                                 build_coastsat_series)
    series = build_coastsat_series(
        [CoastSatDataset(label=f"CoastSat {start}", period_start=start,
                         csv_path=csv_path)], start, module.LOESS_CONFIG)
    cs = series[0]
    match = [w for w in cs["windows"] if w["window"] == module.TARGET_WINDOW]
    if not match:
        raise ValueError(f"window {module.TARGET_WINDOW} not computed")
    loess = dict(zip(np.asarray(match[0]["gis_x"], dtype=int),
                     np.asarray(match[0]["smoothed"], dtype=float)))
    frame = pd.read_csv(csv_path, usecols=["transect_id", "domain_number",
                                           "lrr_m_yr"])
    frame = frame.dropna(subset=["domain_number", "lrr_m_yr"])
    frame["domain"] = frame["domain_number"].astype(int)
    frame["lrr"] = frame["lrr_m_yr"].astype(float)

    # A domain is 500 m holding ~10 transects. Plotted at the domain integer
    # they stack into a vertical column, which reads as one uncertain value
    # rather than as an alongshore gradient. Spread them across the domain in
    # transect order instead. The index is the LAST numeric field of the id --
    # the first is the CoastSat site number and grabbing it sorts by site.
    order = frame["transect_id"].astype(str).str.split("_").str[-1]
    frame["t_index"] = pd.to_numeric(order, errors="coerce")
    frame = frame.sort_values(["domain", "t_index"])
    rank = frame.groupby("domain").cumcount()
    count = frame.groupby("domain")["lrr"].transform("size")
    frame["x"] = frame["domain"] - 0.5 + (rank + 0.5) / count

    sd = frame.groupby("domain")["lrr"].std().to_dict()
    return loess, {int(k): float(v) for k, v in sd.items()}, frame


def run_rates(period, preset, scenario):
    name = f"HAT_{period}_{preset}_{scenario}_groin"
    # Resolved rather than joined by hand -- see the note on the twin of this
    # function in HAT_hindcast_final_figure.py. A hand-built path has no slot
    # for the arm component and reads an arm-scoped run as missing.
    try:
        run_dir = find_run_dir(RAW_RUNS, name, period, preset)
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            f"{exc}  Run HAT_run_all.py first.") from None
    path = run_dir / f"{name}_shoreline_change_rate.csv"
    return pd.read_csv(path).set_index("gis_domain")["lrr_m_yr"], name


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--preset", default="calibBE")
    args = parser.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D

    import dataclasses

    from cascade_pipeline.annotations import (add_geographic_annotations,
                                              annotation_legend_handles)
    from hatteras_site_config import HATTERAS_ANNOTATIONS

    # Label heights are tuned PER FIGURE, not in the shared config, because
    # they depend on this figure's y range and on what is drawn where. The
    # defaults (groin 0.68, piers 0.76) put "Buxton Groin" through the D1-D10
    # transect scatter and "Avon Pier" through both shoreline curves.
    #
    #   Buxton Groin  0.92  above the scatter in (a) and the curves in (b); the
    #                       "Buxton" town label is centred at D7.5 so its box
    #                       clears D5.5 horizontally.
    #   Avon Pier     0.85  the "Avon" town span is ALSO centred on D26, and it
    #                       sits at 0.90, so the pier has to hang below it.
    #   Rodanthe Pier 0.72  the "Rodanthe" village line is at D80 with its label
    #                       at 0.84; dropping the pier clears that, and D79 is
    #                       deeply erosional in both periods so the mid-panel is
    #                       empty there.
    annotations = dataclasses.replace(
        HATTERAS_ANNOTATIONS,
        groin_label_y=0.92,
        piers={"Avon Pier": (26, 0.85), "Rodanthe Pier": (79, 0.72)})

    module = analysis_module()
    reserved = list(module.GROIN_RESERVED_DOMAINS)

    panels = {}
    for key, meta in PERIODS.items():
        csv = (module.P1_COASTSAT_CSV if meta["start"] == 1984
               else module.P2_COASTSAT_CSV)
        loess, sd, transects = loess_and_spread(module, meta["start"], csv)
        model, run_name = run_rates(key, args.preset, meta["scenario"])
        panels[key] = dict(loess=loess, sd=sd, transects=transects,
                           model=model, run=run_name, **meta)

    gis = np.arange(1, 91)
    lo, hi = np.inf, -np.inf
    for d in panels.values():
        for source in (d["loess"], d["model"]):
            values = np.array([source.get(g, np.nan) for g in gis], dtype=float)
            lo = min(lo, np.nanmin(values))
            hi = max(hi, np.nanmax(values))
    pad = 0.10 * (hi - lo)
    # Asymmetric: the top needs room for the place-name strip, the bottom does
    # not, and a symmetric pad on a shared axis wastes a band of panel (b).
    ylim = (lo - pad * 0.45, hi + pad * 1.5)

    plt.rcParams.update({"font.size": 11})
    figure, axes = plt.subplots(2, 1, figsize=(16, 10), sharex=True, sharey=True)

    summary = []
    for axis, (key, d) in zip(axes, panels.items()):
        obs = np.array([d["loess"].get(g, np.nan) for g in gis], dtype=float)
        mod = np.array([d["model"].get(g, np.nan) for g in gis], dtype=float)
        spread = np.array([d["sd"].get(g, np.nan) for g in gis], dtype=float)

        # Communities, village centres, piers, groins and shoal zones, from the
        # shared layer every other cascade_pipeline figure uses -- so this
        # figure cannot disagree with the annotated run plots about where the
        # Buxton groin or the Avon pier is.
        add_geographic_annotations(axis, annotations)

        # The D5-D7 (groin-reserved) and D1/D90 (locked) spans were dropped
        # from this figure: with the geographic annotation layer in place the
        # Buxton area carried four overlapping fills and read as clutter. Both
        # facts still hold and are stated in the caption instead -- and the
        # groin line itself still marks D5.5.
        axis.axvspan(0.5, SKIP_SOUTH + 0.5, facecolor="none",
                     edgecolor="#999999", hatch="\\\\\\", linewidth=0.0,
                     alpha=0.30, zorder=0)

        axis.fill_between(gis, obs - spread, obs + spread, color="#999999",
                          alpha=0.22, zorder=2,
                          label="observed spread (±1 SD of transects)")

        # Individual transects over D1-D10. The LOESS is dashed there because
        # the smoother is unreliable at the end of its range; the scatter is
        # the actual evidence, and it shows the Cape Point spread the smooth
        # curve cannot represent -- roughly 5 m/yr across ten domains.
        south_t = d["transects"][d["transects"]["domain"] <= SKIP_SOUTH]
        axis.plot(south_t["x"], south_t["lrr"], linestyle="none",
                  marker="o", markersize=2.6, color="#4A90D9", alpha=0.62,
                  zorder=5, label="CoastSat transect LRR (D1–D10)")
        axis.fill_between(gis, obs, mod, where=~(np.isnan(obs) | np.isnan(mod)),
                          color=d["colour"], alpha=0.18, zorder=3, label="misfit")

        south = gis <= SKIP_SOUTH
        axis.plot(gis[~south], obs[~south], color="#1A1A1A", linewidth=2.8,
                  zorder=7, label="CoastSat LOESS (10-domain)")
        axis.plot(gis[south], obs[south], color="#1A1A1A", linewidth=2.0,
                  linestyle=(0, (4, 2)), zorder=7,
                  label="LOESS, D1–D10 (excluded by convention)")
        axis.plot(gis, mod, color=d["colour"], linewidth=2.6, zorder=6,
                  label="CASCADE, calibrated")
        axis.axhline(0.0, color="#AAAAAA", linewidth=0.9, zorder=4)

        shared = [g for g in SCORE_DOMAINS
                  if g in d["model"].index and not np.isnan(d["loess"].get(g, np.nan))]
        mm = np.array([d["model"][g] for g in shared])
        oo = np.array([d["loess"][g] for g in shared])
        rmse = float(np.sqrt(np.mean((mm - oo) ** 2)))
        bias = float(np.mean(mm - oo))
        corr = float(np.corrcoef(mm, oo)[0, 1])
        summary.append((d["label"], d["run"], rmse, bias, corr, len(shared)))

        axis.text(0.006, 0.965, f"({d['panel']})", transform=axis.transAxes,
                  fontsize=15, weight="bold", va="top", zorder=12)
        # Right-aligned: left-aligned it ran across D3-D45 at the top and the
        # raised groin and pier labels had nowhere to go.
        axis.text(0.995, 0.965,
                  f"{d['label']}      RMSE {rmse:.3f} m yr⁻¹      "
                  f"bias {bias:+.3f}      r = {corr:.2f}      "
                  f"(D{SKIP_SOUTH + 1}–D89, n = {len(shared)})",
                  transform=axis.transAxes, fontsize=12.5, va="top", ha="right",
                  zorder=12,
                  bbox=dict(facecolor="white", alpha=0.88, edgecolor="none",
                            pad=2.5))

        axis.set_ylabel("Shoreline change rate (m yr⁻¹)\n[+ = seaward]",
                        fontsize=12)
        axis.grid(alpha=0.20, linewidth=0.7)
        axis.set_ylim(ylim)
        axis.tick_params(labelsize=11)

    # End labels from the site config, matching the annotated run plots. The
    # annotation layer already names Buxton, Avon and Tri-Village, so the
    # earlier hand-written place strip only duplicated them.
    axes[0].text(0.0, 1.012, f"← {annotations.low_end_label}",
                 transform=axes[0].transAxes, ha="left", va="bottom",
                 fontsize=10.5, style="italic", color="#444444")
    axes[0].text(1.0, 1.012, f"{annotations.high_end_label} →",
                 transform=axes[0].transAxes, ha="right", va="bottom",
                 fontsize=10.5, style="italic", color="#444444")

    axes[1].set_xlabel("Alongshore domain  (south → north;  500 m per domain)",
                       fontsize=12)
    axes[1].set_xlim(0, 91)
    axes[1].set_xticks([1, 10, 20, 30, 40, 50, 60, 70, 80, 90])

    handles, labels = axes[0].get_legend_handles_labels()
    handles += annotation_legend_handles(annotations)
    # Figure-level, below both panels. Nine entries is too many to sit inside
    # an axis without covering something -- in panel (b) it hid the Avon Shoals
    # label, and in panel (a) every position sat on the CASCADE line.
    figure.legend(handles=handles, loc="lower center", fontsize=9.5, ncol=4,
                  framealpha=0.95, bbox_to_anchor=(0.5, 0.076))

    figure.suptitle("Calibrated CASCADE hindcast against the CoastSat LOESS "
                    "reference, Cape Hatteras", fontsize=15, y=0.983)
    figure.tight_layout(rect=(0, 0.155, 1, 0.963))
    figure.text(
        0.008, 0.008,
        "Configuration: calibrated source/sink field, full management (roadway + beach/dune; nourishment in 2004–2024), groin active at "
        "M = 60 m yr⁻¹, f = 0.6. Observations are the 10-domain LOESS of CoastSat transect rates; the hatched span D1–D10 is drawn dashed because "
        "the project excludes the LOESS there — the smoother is poorly constrained at the end of its range and Cape Point's attachment–detachment "
        "cycle is short-wavelength signal a 5 km smoother removes. Statistics are computed over D11–D89, where the LOESS curve and the calibration "
        "target are identical numbers, so the statistic describes the line that is drawn; the companion figure scores D2–D89 against the spliced "
        "target (raw domain means over D1–D10) and reports 0.526 / 0.558 m yr⁻¹. Source/sink corrections were applied only within a zone set fixed "
        "before calibration; D5–D7 were reserved for the groin, so the misfit there is the groin module's and is not absorbed by the sediment budget.",
        fontsize=7.6, color="#333333", wrap=True)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / "hindcast_calibrated_loess_reference.png"
    figure.savefig(path, dpi=200, facecolor="white")
    plt.close(figure)

    print(f"wrote {path}")
    print(f"  shared y limits: {ylim[0]:.2f} to {ylim[1]:.2f} m/yr")
    for label, run, rmse, bias, corr, n in summary:
        print(f"  {label}  {run}")
        print(f"      RMSE {rmse:.4f}   bias {bias:+.4f}   r {corr:.3f}   n {n}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
