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
PROJECT_BASE_DIR = _HERE.parents[3]
if not (PROJECT_BASE_DIR / "pyproject.toml").exists():
    raise RuntimeError(
        f"CASCADE repo root not found: {PROJECT_BASE_DIR} has no pyproject.toml.")
for _path in (PROJECT_BASE_DIR / "scripts", _HERE.parent):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from hatteras_site_config import HATTERAS_DOMAINS as GEOMETRY  # noqa: E402
from cascade_pipeline.hindcast import implied_interception_m3_yr  # noqa: E402
from HAT_groin_sweep_config import (  # noqa: E402
    GROIN_DOWNDRIFT_GIS,
    GROIN_UPDRIFT_GIS,
)

# The rig pads 11 real domains (D2-D12) with 15 buffer either side, so D5 -> 18
# and D6 -> 19. This is the RIG's convention and differs from production's
# (D5 -> 19, D6 -> 20) -- see HAT_groin_hindcast_1967_2017.py:76.
RIG_BUFFER, RIG_FIRST_GIS, RIG_START_YEAR = 15, 2, 1967

RAW_RUNS = PROJECT_BASE_DIR / "output" / "raw_runs"
RUN = "HAT_1967_2018_edge_calibrated_groin"
DIAGNOSTICS = RAW_RUNS / RUN / f"{RUN}_groin_diagnostics.csv"
SHORELINE = RAW_RUNS / RUN / f"{RUN}_shoreline_matrix.npy"

FIGURE_DIR = PROJECT_BASE_DIR / "output" / "groin_sweep" / "figures"

# Same as HAT_groin_choice_figure.py, so the volumes reconcile: h_b + d_sf.
PROFILE_HEIGHT_M = 1.7 + 22.25
DRIFT_LOW, DRIFT_HIGH = 5.0e5, 7.0e5

# A rig run directory does not name its own parameters -- the sweep writes every
# cell into one run name. Checked against the diagnostics rather than trusted.
EXPECT_M, EXPECT_F = 60.0, 0.6

INK, UP_C, DOWN_C = "#1A1A2E", "#FF8C00", "#2A6FB0"
BAD, MUTED, GRID = "#C2185B", "#777777", "#D5D8DD"


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

    plt.rcParams.update({"font.size": 10, "axes.linewidth": 0.7,
                         "legend.frameon": False, "pdf.fonttype": 42})
    figure, (ax_rate, ax_cum, ax_keep) = plt.subplots(
        3, 1, figsize=(12.0, 11.4), sharex=True,
        gridspec_kw={"height_ratios": [1.0, 1.0, 1.0], "hspace": 0.16})

    millions = FuncFormatter(lambda v, _: f"{v / 1e6:.0f}M")
    thousands = FuncFormatter(lambda v, _: f"{v / 1e3:.0f}k")

    # ---- (a) annual interception against the drift band --------------------
    ax_rate.axhspan(DRIFT_LOW, DRIFT_HIGH, color="#E8F1F8", zorder=0)
    ax_rate.annotate("littoral drift, 5–7 × 10⁵ m³/yr\n(a literature range, not a limit)",
                     xy=(0.985, (DRIFT_LOW + DRIFT_HIGH) / 2),
                     xycoords=("axes fraction", "data"), ha="right",
                     va="center", fontsize=8.6, color="#2A6FB0")
    ax_rate.plot(years, annual_volume, color=UP_C, linewidth=2.4,
                 drawstyle="steps-post", zorder=4)
    ax_rate.axhline(0.0, color="#BBBBBB", linewidth=0.8, zorder=1)
    ax_rate.set_ylabel("Annual transfer across\nthe groin (m³/yr)")
    ax_rate.yaxis.set_major_formatter(thousands)
    ax_rate.set_title(
        "(a)  What the module moves each year — M_eff × 500 m × "
        f"{PROFILE_HEIGHT_M:.1f} m profile",
        fontsize=11, loc="left")
    sound, floor = annual_volume.max(), annual_volume[annual_volume > 0].min()
    ax_rate.annotate(f"M = 60 while sound:  {sound:,.0f} m³/yr",
                     xy=(0.22, 0.55), xycoords="axes fraction",
                     ha="left", va="center", fontsize=9, color=INK)
    ax_rate.annotate(f"after 2003, M·f = 36:  {floor:,.0f} m³/yr",
                     xy=(0.22, 0.38), xycoords="axes fraction",
                     ha="left", va="center", fontsize=9, color=INK)

    # ---- (b) cumulative volume, mirrored ----------------------------------
    ax_cum.fill_between(years, 0, cumulative_volume, color=UP_C, alpha=0.30,
                        zorder=2)
    ax_cum.fill_between(years, 0, -cumulative_volume, color=DOWN_C, alpha=0.30,
                        zorder=2)
    ax_cum.plot(years, cumulative_volume, color=UP_C, linewidth=2.2,
                label="gained by the updrift cell (D6)", zorder=4)
    ax_cum.plot(years, -cumulative_volume, color=DOWN_C, linewidth=2.2,
                label="taken from the downdrift cell (D5)", zorder=4)
    ax_cum.axhline(0.0, color=INK, linewidth=1.2, zorder=3)
    ax_cum.set_ylabel("Cumulative volume\ntransferred (m³)")
    ax_cum.yaxis.set_major_formatter(millions)
    ax_cum.set_title("(b)  The transfer is exactly volume-neutral — "
                     "and the real structure is not",
                     fontsize=11, loc="left")
    ax_cum.legend(loc="upper left", fontsize=9)
    ax_cum.annotate(
        f"±{cumulative_volume[-1] / 1e6:.1f} million m³ by {years[-1]}\n"
        "net across the pair: exactly 0, by construction",
        xy=(0.985, 0.63), xycoords="axes fraction", ha="right", va="center",
        fontsize=8.8, color=INK)
    ax_cum.annotate(
        "THE ASSUMPTION THAT DEPARTS FROM THE FIELD EVIDENCE\n"
        "Observed downdrift extent is 0 m; the model's is 2,500 m. The real\n"
        "structure accretes updrift with no measurable downdrift deficit —\n"
        "the sand comes from the cape, not from D5. Tested: removing the sink\n"
        "halves the fillet and quadruples reach bias, so deleting it is not the fix.",
        xy=(0.015, 0.05), xycoords="axes fraction", ha="left", va="bottom",
        fontsize=8.4, color=BAD)

    # ---- (c) gross against net --------------------------------------------
    ax_keep.plot(years, cumulative_m, color=MUTED, linewidth=2.2,
                 linestyle="--", label="cumulative displacement APPLIED by the module",
                 zorder=3)
    ax_keep.plot(fillet_on_years, fillet_aligned, color=UP_C, linewidth=2.6,
                 label="fillet actually REALISED (D5 − D6)", zorder=4)
    ax_keep.axhline(0.0, color="#BBBBBB", linewidth=0.8, zorder=1)
    ax_keep.set_ylabel("Shoreline displacement (m)")
    ax_keep.set_xlabel("Year")
    ax_keep.set_title("(c)  Almost none of it is retained — "
                      "so M is not a rate of impoundment",
                      fontsize=11, loc="left")
    ax_keep.legend(loc="upper left", fontsize=9)

    retained = 100.0 * fillet_aligned[-1] / cumulative_m[-1]
    ax_keep.annotate(
        f"applied {cumulative_m[-1]:,.0f} m   →   realised "
        f"{fillet_aligned[-1]:.0f} m\n"
        f"retained: {retained:.0f}%   (peak fillet {fillet.max():.0f} m)\n\n"
        "BRIE's alongshore diffusion removes the rest. M is the rate needed to\n"
        "SUSTAIN a fillet against diffusion, not the rate sand is impounded —\n"
        "which is why the 719,000 m³/yr affordability figure is a GROSS\n"
        "restoring rate, not like-for-like against a NET transport budget.",
        xy=(0.985, 0.52), xycoords="axes fraction", ha="right", va="top",
        fontsize=8.8, color=INK)

    for axis in (ax_rate, ax_cum, ax_keep):
        axis.grid(axis="y", color=GRID, linewidth=0.6)
        axis.set_axisbelow(True)
        for side in ("top", "right"):
            axis.spines[side].set_visible(False)

    figure.suptitle(
        f"Sediment budget of the groin module — rig, 1967–{years[-1]}, "
        f"M = {applied_M:g}, f = {applied_f:.1f}",
        fontsize=13, y=0.977)
    figure.text(
        0.012, 0.013,
        "Volumes use the repo's own implied_interception_m3_yr (M × domain "
        "spacing × active profile height), so they reconcile with the "
        "affordability numbers quoted in GROIN_PLAN.md and the run reports. "
        "M is an effective, grid-specific, FIELD-AGGREGATE rate for four "
        "structures inside one 500 m domain against a ~190 m real fillet — not "
        "a sediment flux, and not divisible by four for a per-structure value.",
        fontsize=7.8, color="#555555", wrap=True, va="bottom")

    figure.subplots_adjust(top=0.935, bottom=0.075, left=0.105, right=0.985)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    out = FIGURE_DIR / "sediment_budget.png"
    figure.savefig(out, dpi=200)
    print(f"wrote {out}")
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
