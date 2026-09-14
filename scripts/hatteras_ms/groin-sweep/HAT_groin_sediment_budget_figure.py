#!/usr/bin/env python3
"""The groin module's sediment budget: what it moves, and what it keeps.

WHY THIS FIGURE EXISTS
    `groin_diagnostics.csv` has recorded `cumulative_updrift_m` and
    `cumulative_downdrift_m` every year of every groin run since the module was
    written, and nothing has ever plotted them. Two things that matter for
    reading M were therefore visible nowhere:

      1. THE DIPOLE IS VOLUME-NEUTRAL AND THE REAL STRUCTURE IS NOT. The module
         takes exactly as much from the downdrift cell as it gives the updrift
         one. The real Buxton structure accretes updrift with NO measurable
         downdrift deficit -- observed downdrift extent 0 m against the model's
         2,500 m -- because the sand comes from the cape, not from D5. This is
         the module's largest structural assumption and it had no figure.

      2. ALMOST NONE OF WHAT THE DIPOLE INJECTS IS RETAINED. Over the rig's 50
         years at M = 60, f = 0.6 the module applies ~2,400 m of cumulative
         one-sided displacement and holds a fillet of ~69 m. BRIE's alongshore
         diffusion removes the rest. So M is NOT the rate at which sand is
         impounded -- it is the rate needed to SUSTAIN a fillet against
         diffusion, which is a much larger number.

    Panel (c) is the reason this figure is worth having. It reframes the
    affordability comparison that GROIN_PLAN.md and the run reports both make:
    719,000 m3/yr at M = 60 against a 5-7e5 m3/yr littoral drift is a GROSS
    restoring rate set against a NET transport budget, and they are not like
    for like. That does not make the comparison wrong -- it is a deliberate,
    documented diagnostic -- but it does mean "marginally above the drift band"
    should not be read as "impounds more sand than the coast carries."

WHAT IS PLOTTED
    (a) ANNUAL INTERCEPTION. M_eff each year converted to a volume by the
        repo's own `implied_interception_m3_yr` (M * dy * profile height),
        against the shaded 5-7e5 m3/yr littoral drift band. Shows the
        deterioration schedule carrying the module from marginally above the
        band down into it.
    (b) CUMULATIVE VOLUME, updrift and downdrift as exact mirror images. The
        symmetry IS the assumption; the annotation is where it departs from the
        field evidence.
    (c) GROSS AGAINST NET. Cumulative applied displacement against the fillet
        actually realised, same axis, same units.

    Volumes use the repo's conversion and the same profile height as
    HAT_groin_choice_figure.py, so the numbers reconcile with the affordability
    figures quoted elsewhere.

Usage:
    python HAT_groin_sediment_budget_figure.py

Writes output/groin_sweep/figures/sediment_budget.png

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
PROJECT_BASE_DIR = next(_p for _p in _HERE.parents
                        if (_p / "pyproject.toml").exists())
if not (PROJECT_BASE_DIR / "pyproject.toml").exists():
    raise RuntimeError(
        f"CASCADE repo root not found: {PROJECT_BASE_DIR} has no pyproject.toml.")
for _path in (PROJECT_BASE_DIR / "scripts", _HERE.parent):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from hatteras_site_config import HATTERAS_DOMAINS as GEOMETRY  # noqa: E402
from cascade_pipeline.hindcast import implied_interception_m3_yr  # noqa: E402
from cascade_pipeline.run_layout import resolve  # noqa: E402
from hat_figure_style import (apply_style, C, INK, INK_MUTED,  # noqa: E402
                              caption, figsize, open_frame, save, _title)
from HAT_groin_sweep_config import (  # noqa: E402
    GROIN_DOWNDRIFT_GIS,
    GROIN_UPDRIFT_GIS,
)

# The rig pads 11 real domains (D2-D12) with 15 buffer either side, so D5 -> 18
# and D6 -> 19. This is the RIG's convention and differs from production's
# (D5 -> 19, D6 -> 20) -- see HAT_groin_hindcast_1967_2017.py:76.
RIG_BUFFER, RIG_FIRST_GIS, RIG_START_YEAR = 15, 2, 1967

# The rig lives in output/rig_runs/, not output/raw_runs/ (moved 2026-08-31).
# It is a DIFFERENT GRID -- 41 domains against production's 120 -- and M is
# grid-specific, so mixing the two invited quoting a rig number as a
# production one. raw_runs is production only, and run_index.csv covers it.
RAW_RUNS = PROJECT_BASE_DIR / "output" / "rig_runs"
RUN = "HAT_1967_2018_edge_calibrated_groin"
DIAGNOSTICS = resolve(RAW_RUNS / RUN, "groin_csv", RUN)
SHORELINE = resolve(RAW_RUNS / RUN, "matrix", RUN)

FIGURE_DIR = PROJECT_BASE_DIR / "output" / "groin_sweep" / "figures"

# Same as HAT_groin_choice_figure.py, so the volumes reconcile: h_b + d_sf.
PROFILE_HEIGHT_M = 1.7 + 22.25
DRIFT_LOW, DRIFT_HIGH = 5.0e5, 7.0e5

# A rig run directory does not name its own parameters -- the sweep writes every
# cell into one run name. Checked against the diagnostics rather than trusted.
EXPECT_M, EXPECT_F = 60.0, 0.6

# The two sides of the structure, in the house semantics: the sheltered
# updrift cell is the ACCENT (the thing the module does), the downdrift
# cell it takes from is BASE grey. INK/UP_C/DOWN_C/BAD/MUTED/GRID were a
# navy, an orange, a blue, a pink and two greys chosen in this file.
UP_C, DOWN_C = C["ACCENT"], C["BASE"]


def _rig_pad(gis_id: int) -> int:
    return RIG_BUFFER + (gis_id - RIG_FIRST_GIS)


def main() -> None:
    for path in (DIAGNOSTICS, SHORELINE):
        if not path.exists():
            raise SystemExit(
                f"missing rig output: {path}\n"
                "Run HAT_groin_hindcast_1967_2017.py first (run key 'groin').")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter

    diagnostics = pd.read_csv(DIAGNOSTICS)
    active = diagnostics[diagnostics["groin_active"]]
    applied_M = float(active["trapping_rate_applied_m_yr"].max())
    applied_f = float(active["trapping_rate_applied_m_yr"].min()) / applied_M
    if not (np.isclose(applied_M, EXPECT_M)
            and np.isclose(applied_f, EXPECT_F, atol=0.02)):
        raise SystemExit(
            f"{DIAGNOSTICS.name} reports M = {applied_M:g}, f = {applied_f:.2f}, "
            f"not the expected {EXPECT_M:g} / {EXPECT_F:g}. This run directory "
            "holds a different cell. Re-run HAT_groin_hindcast_1967_2017.py.")

    years = diagnostics["model_year"].to_numpy()
    rate = diagnostics["trapping_rate_applied_m_yr"].to_numpy()
    annual_volume = np.array(
        [implied_interception_m3_yr(r, PROFILE_HEIGHT_M, GEOMETRY) for r in rate])
    # Sign convention in the CSV: updrift is negative (seaward in the model's
    # landward-positive frame), downdrift positive. Magnitudes are identical by
    # construction -- that identity is the point of panel (b).
    cumulative_m = diagnostics["cumulative_downdrift_m"].to_numpy()
    cumulative_volume = cumulative_m * GEOMETRY.domain_spacing_m * PROFILE_HEIGHT_M

    matrix = np.load(SHORELINE)
    fillet = (matrix[:, _rig_pad(GROIN_DOWNDRIFT_GIS)]
              - matrix[:, _rig_pad(GROIN_UPDRIFT_GIS)])
    fillet = fillet - fillet[0]
    fillet_years = RIG_START_YEAR + np.arange(matrix.shape[0])
    # The diagnostics stop one year before the shoreline matrix; align on year.
    keep = np.isin(fillet_years, years)
    fillet_on_years, fillet_aligned = fillet_years[keep], fillet[keep]

    apply_style()
    figure, (ax_rate, ax_cum, ax_keep) = plt.subplots(
        3, 1, figsize=figsize("double", aspect=1.05), sharex=True,
        gridspec_kw={"height_ratios": [1.0, 1.0, 1.0]},
        constrained_layout=True)

    millions = FuncFormatter(lambda v, _: f"{v / 1e6:.0f}M")
    thousands = FuncFormatter(lambda v, _: f"{v / 1e3:.0f}k")

    # ---- (a) annual interception against the drift band --------------------
    ax_rate.axhspan(DRIFT_LOW, DRIFT_HIGH, color="0.94", zorder=0)
    ax_rate.annotate("littoral drift, 5 to 7 × 10⁵ m³/yr",
                     xy=(0.985, (DRIFT_LOW + DRIFT_HIGH) / 2),
                     xycoords=("axes fraction", "data"), ha="right",
                     va="center", fontsize=7, color=INK_MUTED)
    ax_rate.plot(years, annual_volume, color=UP_C, linewidth=1.6,
                 drawstyle="steps-post", zorder=4)
    ax_rate.axhline(0.0, color=INK_MUTED, linewidth=0.8, zorder=1)
    ax_rate.set_ylabel("annual transfer across\nthe groin (m³/yr)")
    ax_rate.yaxis.set_major_formatter(thousands)
    _title(ax_rate, 0, "what the module moves each year")
    sound, floor = annual_volume.max(), annual_volume[annual_volume > 0].min()
    ax_rate.annotate(f"while sound, {sound:,.0f} m³/yr",
                     xy=(0.22, 0.58), xycoords="axes fraction",
                     ha="left", va="center", fontsize=7.5, color=INK)
    ax_rate.annotate(f"after 2003, {floor:,.0f} m³/yr",
                     xy=(0.22, 0.42), xycoords="axes fraction",
                     ha="left", va="center", fontsize=7.5, color=INK)

    # ---- (b) cumulative volume, mirrored ----------------------------------
    ax_cum.fill_between(years, 0, cumulative_volume, color=UP_C, alpha=0.25,
                        linewidth=0, zorder=2)
    ax_cum.fill_between(years, 0, -cumulative_volume, color=DOWN_C, alpha=0.25,
                        linewidth=0, zorder=2)
    ax_cum.plot(years, cumulative_volume, color=UP_C, linewidth=1.6,
                label="gained by the updrift cell, D6", zorder=4)
    ax_cum.plot(years, -cumulative_volume, color=DOWN_C, linewidth=1.6,
                label="taken from the downdrift cell, D5", zorder=4)
    ax_cum.axhline(0.0, color=INK, linewidth=0.8, zorder=3)
    ax_cum.set_ylabel("cumulative volume\ntransferred (m³)")
    ax_cum.yaxis.set_major_formatter(millions)
    _title(ax_cum, 1, "the transfer is exactly volume-neutral")
    ax_cum.legend(loc="upper left", fontsize=7.5)
    ax_cum.annotate(
        f"±{cumulative_volume[-1] / 1e6:.1f} million m³ by {years[-1]};"
        " net across the pair, exactly zero",
        xy=(0.985, 0.62), xycoords="axes fraction", ha="right", va="center",
        fontsize=7.5, color=INK)

    # ---- (c) gross against net --------------------------------------------
    ax_keep.plot(years, cumulative_m, color=C["BASE"], linewidth=1.4,
                 linestyle="--",
                 label="cumulative displacement applied by the module",
                 zorder=3)
    ax_keep.plot(fillet_on_years, fillet_aligned, color=UP_C, linewidth=1.8,
                 label="fillet actually realised, D5 − D6", zorder=4)
    ax_keep.axhline(0.0, color=INK_MUTED, linewidth=0.8, zorder=1)
    ax_keep.set_ylabel("shoreline displacement (m)")
    ax_keep.set_xlabel("year")
    _title(ax_keep, 2, "almost none of it is retained")
    ax_keep.legend(loc="upper left", fontsize=7.5)

    retained = 100.0 * fillet_aligned[-1] / cumulative_m[-1]
    ax_keep.annotate(
        f"applied {cumulative_m[-1]:,.0f} m, realised "
        f"{fillet_aligned[-1]:.0f} m: {retained:.0f}% retained",
        xy=(0.985, 0.52), xycoords="axes fraction", ha="right", va="top",
        fontsize=7.5, color=INK)

    for axis in (ax_rate, ax_cum, ax_keep):
        axis.grid(axis="y")
        axis.set_axisbelow(True)
        open_frame(axis)

    caption(figure,
            "The sediment budget of the groin module on the 1967 rig, "
            "1967 to {end}, at M = {M:g} and f = {f:.1f}. (a) What the module "
            "moves each year: the effective trapping rate times the 500 m "
            "domain spacing times a {h:.1f} m active profile, which is "
            "{sound:,.0f} m³/yr while the structure is sound and {floor:,.0f} "
            "m³/yr after the 2003 storm. The band is a literature range for "
            "the littoral drift, 5 to 7 × 10⁵ m³/yr — a comparison, not a "
            "limit the module enforces. (b) The transfer is exactly "
            "volume-neutral: the updrift cell gains what the downdrift cell "
            "loses, ±{cum:.1f} million m³ by {end}, netting to zero by "
            "construction. THAT IS THE ASSUMPTION THAT DEPARTS FROM THE FIELD "
            "EVIDENCE. The observed downdrift extent is 0 m against the "
            "model's 2,500 m: the real structure accretes updrift with no "
            "measurable downdrift deficit, because the sand comes from the "
            "cape and not from D5. Removing the sink was tested and is not the "
            "fix — it halves the fillet and quadruples the reach bias. "
            "(c) Almost none of what is applied is retained: {applied:,.0f} m "
            "of displacement produces {real:.0f} m of fillet, {ret:.0f}%, "
            "because BRIE's alongshore diffusion removes the rest. So M is the "
            "rate needed to SUSTAIN a fillet against diffusion, not a rate at "
            "which sand is impounded, and the affordability figure it implies "
            "is a GROSS restoring rate rather than like-for-like against a net "
            "transport budget. Volumes come from the repo's own "
            "implied_interception_m3_yr, so they reconcile with the numbers in "
            "GROIN_PLAN.md and the run reports. M is an effective, "
            "grid-specific, field-aggregate rate for four structures inside "
            "one 500 m domain against a real fillet of about 190 m: not a "
            "sediment flux, and not divisible by four for a per-structure "
            "value."
            .format(end=years[-1], M=applied_M, f=applied_f,
                    h=PROFILE_HEIGHT_M, sound=sound, floor=floor,
                    cum=cumulative_volume[-1] / 1e6,
                    applied=cumulative_m[-1], real=fillet_aligned[-1],
                    ret=retained))

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    written = save(figure, FIGURE_DIR / "sediment_budget.png", close=True)
    print(f"wrote {written[0]}")
    print(f"  annual while sound   {sound:,.0f} m3/yr "
          f"({100 * sound / DRIFT_HIGH:.0f}% of the upper drift bound)")
    print(f"  annual after 2003    {floor:,.0f} m3/yr "
          f"({100 * floor / DRIFT_HIGH:.0f}%)")
    print(f"  cumulative by {years[-1]}   {cumulative_volume[-1]:,.0f} m3 "
          f"(±, net 0)")
    print(f"  applied {cumulative_m[-1]:,.0f} m -> realised "
          f"{fillet_aligned[-1]:.1f} m = {retained:.1f}% retained")


if __name__ == "__main__":
    main()
