"""
Net change 2009-2025 with the 2017 Buxton fill shifted two domains north, in the run-figure layout.

    python scripts/hatteras_ms/experiments/HAT_buxton_2017_shift_north_plot.py

Draws the shifted run the way the runner draws every run's net-change figure
(cascade_pipeline.plotting.net_change_comparison: same target, smoothing, annotation layer
and fixed ±100 m axis), then adds the reported-footprint run as a dashed line and both
2017 footprints as bars. Reads the two runs written by HAT_buxton_2017_shift_north.py and
its tables/scores.csv for the legend RMSE.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-08
"""
from __future__ import annotations

import dataclasses
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

_HERE = Path(__file__).resolve()
PROJECT_ROOT = next(_p for _p in _HERE.parents if (_p / "pyproject.toml").exists())
sys.path[:0] = [str(PROJECT_ROOT / "scripts"), str(_HERE.parent)]
import HAT_buxton_2017_shift_north as X  # noqa: E402
from cascade_pipeline.coastsat_lowess import CoastSatDataset, build_coastsat_series  # noqa: E402
from cascade_pipeline.plotting import net_change_comparison as N  # noqa: E402
from cascade_pipeline.plotting.rate_comparison import (DEFAULT_RATE_COMPARISON, _publication_legend,  # noqa: E402
                                                       annotation_legend_handles, plot_annotated_rate_comparison)
from cascade_pipeline.run_info import RunInfo  # noqa: E402
from cascade_pipeline.shoreline import compute_change_rate  # noqa: E402
from site_layer.hat_figure_style import record_caption  # noqa: E402
from site_layer.hat_observed_rates import NET_CHANGE_CENTRES, NET_CHANGE_WINDOWS, WINDOW_ROLE, net_change_dir  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS as A, HATTERAS_DOMAINS as D  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
OUT = X.EXP_DIR / "figures" / "buxton_2017_reported_vs_shifted_north_net_change_2009_2025.png"
REPORTED_2017 = (6, 16)
REPORTED_DASH = (0, (3, 1.6))
FOOTPRINT_Y = {"reported": 93.0, "shifted": 87.0}   # bar heights, m, above the town labels
FOOTPRINT_COLOR = {"reported": "#8C8C8C", "shifted": A.model_color}
# -----------------------------------------------------------------------------


def run_info(rd):
    meta = json.loads(next(Path(rd).glob("*_run_metadata.json")).read_text(encoding="utf-8"))
    p = meta["period"]
    return RunInfo(run_name=Path(rd).name, run_dir=str(rd), start_year=int(p["start_year"]),
                   end_year=int(p["end_year"]), Hs=float(meta["wave climate"]["wave_height_m"]),
                   flip_sign_model=True, background_erosion_on=True,
                   wave_climate=meta.get("scenario", {}).get("wave climate"))


def position(rd):
    m = np.load(next(Path(rd).glob("*_shoreline_matrix.npy")))
    return m, compute_change_rate(m, span_years=1, flip_sign=True)


def main():
    rd = {v: X.run_dir(v) for v in X.VARIANTS}
    run = run_info(rd["shifted"])
    window = (run.start_year, run.end_year)
    m_shift, pos = position(rd["shifted"])
    _, pos_rep = position(rd["reported"])
    sc = pd.read_csv(X.EXP_DIR / "tables" / "scores.csv").set_index(["variant", "area"])
    interior = lambda v, k: sc.loc[(v, "interior GIS 2-89, smoothed"), k]  # noqa: E731
    shifted = (REPORTED_2017[0] + X.SHIFT, REPORTED_2017[1] + X.SHIFT)

    # The run figure, exactly as net_change_comparison draws it
    tag = f"{window[0]}_{window[1]}"
    dataset = CoastSatDataset(
        label=f"CoastSat net change ({window[0]}–{window[1]})", period_start=window[0],
        csv_path=str(net_change_dir(*window) / f"transect_net_change_{tag}.csv"), rate_col="net_change_m")
    cs_series = build_coastsat_series([dataset], active_period_start=window[0],
                                      lowess_config=N.LOWESS_CONFIG, domains=D)
    source_sink, groin, _ = N.setup_text(run.run_name)
    role = WINDOW_ROLE.get(window, "").replace(" (held out)", "")
    config = dataclasses.replace(
        DEFAULT_RATE_COMPARISON, quantity="position", publication_text=True,
        show_features=True, plot_domain_means=True, raw_lrr_southern_only=True,
        observed_label="CoastSat (7-domain LOWESS)",
        model_label=f"CASCADE, 2017 fill on GIS {shifted[0]}–{shifted[1]} (RMSE {interior('shifted', 'rmse_m'):.1f} m)",
        title=f"Modelled and observed shoreline position change, {window[0]}–{window[1]}",
        ylim=(-N.Y_HALF, N.Y_HALF),
        subtitle=" · ".join((role, source_sink, groin,
                             f"experiment: 2017 Buxton fill moved {X.SHIFT} domains north")))
    smoothed = N.smooth_model(pos)
    fig, ax = plot_annotated_rate_comparison(
        smoothed, cs_series, run, estimator=None, domains=D, annotations=A,
        lowess_config=N.LOWESS_CONFIG, config=config, save_path=None, show=False, change_rate_raw=pos)

    # The reported-footprint run, smoothed the same way, dashed
    gis = np.arange(D.first_gis_id, D.last_gis_id + 1)
    real = slice(D.start_real_index, D.end_real_index)
    ax.plot(gis, N.smooth_model(pos_rep)[real], color=A.model_color, lw=1.8, ls=REPORTED_DASH,
            alpha=0.85, zorder=9)

    # The two 2017 footprints
    for v, (lo, hi) in (("reported", REPORTED_2017), ("shifted", shifted)):
        ax.plot([lo - 0.5, hi + 0.5], [FOOTPRINT_Y[v]] * 2, color=FOOTPRINT_COLOR[v], lw=3.2,
                solid_capstyle="butt", zorder=12)
        ax.text(hi + 0.9, FOOTPRINT_Y[v], f"2017 fill, {v} (GIS {lo}–{hi})", fontsize=7,
                color=FOOTPRINT_COLOR[v], va="center", ha="left", zorder=12)

    offaxis = N._mark_offaxis(ax, smoothed, pos, cs_series)
    for leg in list(fig.legends):
        leg.remove()
    extra = [Line2D([0], [0], color=A.model_color, lw=1.8, ls=REPORTED_DASH, alpha=0.85,
                    label=f"CASCADE, 2017 fill on GIS {REPORTED_2017[0]}–{REPORTED_2017[1]} as reported "
                          f"(RMSE {interior('reported', 'rmse_m'):.1f} m)")]
    _publication_legend(fig, config, A, N.LOWESS_CONFIG, extra=extra + list(annotation_legend_handles(A)),
                        ncol=4, extra_model_raw=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=300, bbox_inches="tight", facecolor="white")
    (OUT.parent / "supporting").mkdir(exist_ok=True)
    fig.savefig(OUT.parent / "supporting" / OUT.with_suffix(".pdf").name, bbox_inches="tight", facecolor="white")
    plt.close("all")

    (s0, s1), (e0, e1) = NET_CHANGE_WINDOWS[window]
    c0, c1 = NET_CHANGE_CENTRES[window]
    years = m_shift.shape[0] - 1
    end_note = (f" CoastSat ends on {N.COASTSAT_END}, so that mean covers {e0} to {N.COASTSAT_END}."
                if e1 > N.COASTSAT_END else "")
    record_caption(OUT, (
        f"Experiment: the 2017 Buxton fill moved {X.SHIFT} domains north. Shoreline position change along "
        f"Hatteras Island, {window[0]}–{window[1]} ({role.lower()}), from GIS domain 1 (Cape Point) to 90 "
        f"(Pea Island); positive is seaward. Both model runs are the {years}-year test run (1 Jan "
        f"{window[0] + years} minus the start) with {source_sink}, {groin}, road maintenance, beach and dune "
        f"management and nourishment, relocations off. Solid orange: the 2017 fill on GIS {shifted[0]}–{shifted[1]}; "
        f"thick smoothed with a 7-domain LOWESS north of domain 10 (unsmoothed south of it and at domain 90), "
        f"thin each domain's unsmoothed value. Dashed orange: the same run with the 2017 fill on its reported "
        f"footprint, GIS {REPORTED_2017[0]}–{REPORTED_2017[1]}, smoothed the same way; it is the matrix run "
        f"exactly. The fill's volume (2.6 M cy) and width are unchanged, so the sand per metre is the same; the "
        f"2022 Buxton fill stays on GIS 6–16 in both. Bars near the top mark the two 2017 footprints. Legend "
        f"RMSE is the interior score, GIS 2–89, both sides smoothed; bias is "
        f"{interior('reported', 'bias_m'):+.1f} m (reported) and {interior('shifted', 'bias_m'):+.1f} m "
        f"(shifted), r {interior('reported', 'r'):.2f} and {interior('shifted', 'r'):.2f}. Outside GIS 3–20 "
        f"the runs are identical and the solid line hides the dashed one. "
        f"Observed (CoastSat, blue): per transect, the mean shoreline position over {e0} to {e1} (centred on "
        f"{c1}) minus the mean over {s0} to {s1} (centred on {c0}).{end_note} Thin line: the mean of each "
        f"500 m domain's transects. Thick line: the transect values smoothed with a 7-domain LOWESS north of "
        f"domain 10; south of it the individual transects are shown as dots. Blue shading: the smoothed "
        f"observed change against zero. Shaded bands: communities; dashed lines: village centres; hatched: "
        f"shoal zones; dash-dot: piers; dotted: Buxton groin. The y axis is fixed at ±{N.Y_HALF:g} m as on "
        f"every run figure." + offaxis + (f" Wave climate: {run.wave_climate}." if run.wave_climate else "")))
    print(OUT)


if __name__ == "__main__":
    main()
