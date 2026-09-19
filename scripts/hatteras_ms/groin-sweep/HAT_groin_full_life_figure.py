#!/usr/bin/env python3
"""The groin module over the STRUCTURE'S WHOLE LIFE, against every survey.

WHY THIS FIGURE EXISTS
    `HAT_groin_timeseries_check.py` plots the two hindcast windows, and both
    of them start 15 years after the groin went in. They therefore show the
    fillet's DECAY and never its CREATION -- which is why they cannot see f,
    and why they make the module look worse than it is. The 1967-2018 rig is
    the only window that contains the build phase, the 1996 repair, the 2003
    storm damage and the decline that follows. This figure is that window.

    It is also the figure that justifies f = 0.6: the deterioration ramp is
    visible here and nowhere else.

WHAT IS PLOTTED
    (a) FILLET THROUGH TIME. The surveyed fillet (D5 - D6 offset against the
        fixed 1967 datum, from 24 dated wet/dry surveys) as markers -- markers
        only, because the record is irregular and a joining line would imply
        samples that do not exist -- against the rig's modelled fillet. The
        structure's documented timeline is marked on top.

    (b) WHAT THE MODULE WAS DOING. The trapping rate the module actually
        applied each year, read from the run's own groin_diagnostics.csv
        rather than recomputed. This is the schedule f parameterises: zero
        before install, M while the structure is sound, a linear ramp down
        from the 1996 repair to the 2003 damage, then a hold at M*f.

    Reading the two together is the point. Panel (b) explains the shape of
    the model curve in (a), and the observed peak in (a) at 2004 is what
    fixes the end of the ramp in (b).

HOW TO READ THE RESIDUAL
    The module reproduces the SIGN and the TIMING of the build, and it
    undershoots the AMPLITUDE. That is expected and documented: no admissible
    M matches the fillet on this grid, because the real fillet is ~190 m wide
    against a 500 m domain and the dipole is volume-neutral where the real
    structure is not (observed downdrift extent 0 m, the model's 2,500 m).
    The gap is the part the source/sink calibration and the Cape Point
    dynamics absorb -- not a failed fit.

    Note also that the rig runs a 1967 window off 1984 topography
    (RIG_TOPO_PRODUCT = "1984-start"), a deliberate anachronism accepted
    because the target is a shoreline OFFSET rather than an elevation.

Usage:
    python HAT_groin_full_life_figure.py

Writes output/calibration/groin/figures/full_life_1967_2017.png

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import re
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

from cascade_pipeline.run_layout import resolve  # noqa: E402
from site_layer.hat_figure_style import (apply_style, C, INK, INK_MUTED,  # noqa: E402
                              caption, figsize, open_frame, save, _title)
from HAT_groin_sweep_config import (  # noqa: E402
    GROIN_SWEEP_ROOT,
    GROIN_DOWNDRIFT_GIS,
    GROIN_UPDRIFT_GIS,
    WETDRY_CHANGE_TABLE,
)

# ---------------------------------------------------------------------------
# The rig, and its own index convention.
# ---------------------------------------------------------------------------
# The rig pads 11 real domains (D2-D12) with 15 buffer domains either side, so
# D2 -> 15 and D5 -> 18, D6 -> 19. This is _gis_to_pad() in
# HAT_groin_hindcast_1967_2017.py:76, restated rather than imported because
# importing that module builds a CASCADE run.
RIG_BUFFER = 15
RIG_FIRST_GIS = 2
RIG_START_YEAR = 1967

# The rig lives in output/calibration/groin_rig/, not output/raw_runs/ (moved 2026-08-31;
# under calibration/ since 2026-09-18).
# It is a DIFFERENT GRID -- 41 domains against production's 120 -- and M is
# grid-specific, so mixing the two invited quoting a rig number as a
# production one. raw_runs is production only, and run_index.csv covers it.
RAW_RUNS = PROJECT_BASE_DIR / "output" / "calibration" / "groin_rig"
GROIN_RUN = "HAT_1967_2018_edge_calibrated_groin"
NO_GROIN_RUN = "HAT_1967_2018_edge_calibrated_no_groin"


def _matrix(run_name):
    return resolve(RAW_RUNS / run_name, "matrix", run_name)


SHORELINE = _matrix(GROIN_RUN)
NO_GROIN_SHORELINE = _matrix(NO_GROIN_RUN)
DIAGNOSTICS = resolve(RAW_RUNS / GROIN_RUN, "groin_csv", GROIN_RUN)

# A rig run directory does NOT name its own parameters -- the sweep writes every
# cell into one run name, so whatever survives is the last cell that finished.
# On 2026-08-30 this directory was found holding an UNSTABLE M = 70 cell whose
# fillet ran away to 444 m, while being named as though it were the calibrated
# run. So the applied rate is read back from the diagnostics and checked here
# rather than trusted from the label.
EXPECT_M, EXPECT_F = 60.0, 0.6

FIGURE_DIR = GROIN_SWEEP_ROOT / "figures"

# House colours (2026-09-11): the surveys are INK, the run under test the
# ACCENT, the groin-off run BASE grey, and the structure's dated events are
# guide lines in muted ink rather than a fourth hue. These were near-black,
# an orange, a blue, a pink and a grey chosen in this file.
OBSERVED_COLOR, MODEL_COLOR, RATE_COLOR = INK, C["ACCENT"], C["ACCENT"]
EVENT_COLOR = INK_MUTED

# Documented structure history -- GROIN_PLAN.md section 1.
EVENTS = [
    (1969, "installed"),
    (1996, "last repaired"),
    (2003, "storm damage"),
    (2004, "fillet peaks"),
]


def _rig_pad(gis_id: int) -> int:
    """D2->15, D5->18, D6->19, D12->25."""
    return RIG_BUFFER + (gis_id - RIG_FIRST_GIS)


def observed_fillet_by_year() -> dict:
    """Surveyed fillet against the fixed 1967 datum, {year: metres}.

    Same convention as HAT_groin_timeseries_check.py: downdrift minus
    updrift, so a rising curve means the updrift side is holding while the
    downdrift side retreats -- what a groin builds.
    """
    frame = pd.read_csv(WETDRY_CHANGE_TABLE).set_index("Domain_ID")
    out: dict = {}
    for column in frame.columns:
        # The SECOND year is the survey year; the first is the 1967 datum.
        match = re.match(r"change_from_wetdry_1967_wetdry_(\d{4})", column)
        if not match:
            continue
        up = frame.loc[GROIN_UPDRIFT_GIS, column]
        down = frame.loc[GROIN_DOWNDRIFT_GIS, column]
        if pd.isna(up) or pd.isna(down):
            continue
        out.setdefault(int(match.group(1)), []).append(float(down - up))
    return {year: float(np.mean(v)) for year, v in sorted(out.items())}


def main() -> None:
    for path in (SHORELINE, DIAGNOSTICS):
        if not path.exists():
            raise SystemExit(
                f"missing rig output: {path}\n"
                "Run HAT_groin_hindcast_1967_2017.py first "
                "(run key 'groin').")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def fillet(path):
        matrix = np.load(path)
        offset = (matrix[:, _rig_pad(GROIN_DOWNDRIFT_GIS)]
                  - matrix[:, _rig_pad(GROIN_UPDRIFT_GIS)])
        return (RIG_START_YEAR + np.arange(matrix.shape[0]),
                offset - offset[0])

    years, modelled = fillet(SHORELINE)
    no_groin = (fillet(NO_GROIN_SHORELINE)[1]
                if NO_GROIN_SHORELINE.exists() else None)

    observed = observed_fillet_by_year()
    in_window = {y: v for y, v in observed.items()
                 if years[0] <= y <= years[-1]}
    obs_years = np.array(sorted(in_window))
    obs_vals = np.array([in_window[y] for y in obs_years])

    diagnostics = pd.read_csv(DIAGNOSTICS)
    active = diagnostics[diagnostics["groin_active"]]
    applied_M = float(active["trapping_rate_applied_m_yr"].max())
    applied_f = float(active["trapping_rate_applied_m_yr"].min()) / applied_M
    if not (np.isclose(applied_M, EXPECT_M)
            and np.isclose(applied_f, EXPECT_F, atol=0.02)):
        raise SystemExit(
            f"{DIAGNOSTICS.name} reports M = {applied_M:g}, f = {applied_f:.2f}, "
            f"not the expected {EXPECT_M:g} / {EXPECT_F:g}. This run directory "
            "holds a different cell -- see the note at the head of this file. "
            "Re-run HAT_groin_hindcast_1967_2017.py before plotting.")

    apply_style()
    figure, (ax_fillet, ax_rate) = plt.subplots(
        2, 1, figsize=figsize("double", aspect=0.66), sharex=True,
        gridspec_kw={"height_ratios": [2.4, 1.0]}, constrained_layout=True)

    # ---- (a) fillet ------------------------------------------------------
    for year, label in EVENTS:
        ax_fillet.axvline(year, color=EVENT_COLOR, linewidth=0.8,
                          linestyle=(0, (1, 2)), zorder=2)
        # Inside the axes: above the top edge these ran through the panel
        # title at the printed width.
        ax_fillet.annotate(label, xy=(year, 0.985),
                           xycoords=("data", "axes fraction"),
                           rotation=90, ha="right", va="top",
                           fontsize=7.0, color=EVENT_COLOR)

    ax_fillet.axhline(0.0, color=INK_MUTED, linewidth=0.8,
                      linestyle=(0, (4, 3)), zorder=1)
    if no_groin is not None:
        ax_fillet.plot(years, no_groin, color=C["BASE"], linestyle=":",
                       linewidth=1.4, label="modelled, groin off", zorder=3)
    ax_fillet.plot(years, modelled, color=MODEL_COLOR, linewidth=1.8,
                   label=f"modelled, M {applied_M:g}, f {applied_f:.1f}",
                   zorder=4)
    ax_fillet.plot(obs_years, obs_vals, marker="o", markersize=3.6,
                   linestyle="none", color=OBSERVED_COLOR,
                   label=f"surveyed fillet, {len(obs_years)} dates", zorder=5)

    peak_year = int(obs_years[int(np.argmax(obs_vals))])
    ax_fillet.set_ylabel("fillet, D5 − D6 offset since 1967 (m)\n"
                         "positive is the updrift side holding")
    _title(ax_fillet, 0, "the modelled fillet against every survey")
    ax_fillet.grid(axis="y")
    ax_fillet.set_axisbelow(True)
    ax_fillet.legend(loc="upper left", fontsize=7.5)
    open_frame(ax_fillet)

    final_obs = float(obs_vals[-1])
    final_mod = float(modelled[-1])
    ax_fillet.annotate(
        f"observed peak {obs_vals.max():.0f} m in {peak_year};"
        f" modelled peak {modelled.max():.0f} m",
        xy=(0.985, 0.08), xycoords="axes fraction", ha="right", va="bottom",
        fontsize=7.5, color=INK_MUTED)

    # ---- (b) what the module applied -------------------------------------
    ax_rate.plot(diagnostics["model_year"],
                 diagnostics["trapping_rate_m_yr_applied"]
                 if "trapping_rate_m_yr_applied" in diagnostics
                 else diagnostics["trapping_rate_applied_m_yr"],
                 color=RATE_COLOR, linewidth=1.6, drawstyle="steps-post")
    for year, _ in EVENTS:
        ax_rate.axvline(year, color=EVENT_COLOR, linewidth=0.8,
                        linestyle=(0, (1, 2)), zorder=2)
    ax_rate.set_ylabel("applied trapping\nrate (m/yr)")
    ax_rate.set_xlabel("year")
    ax_rate.grid(axis="y")
    ax_rate.set_axisbelow(True)
    open_frame(ax_rate)
    _title(ax_rate, 1, "the rate the module applied")

    caption(figure,
            "The Buxton groin over its whole life on the 1967 rig, 1967 to "
            "2017, at M = {M:g} and f = {f:.1f}. (a) The modelled fillet "
            "against every survey in the window, with the groin-off run for "
            "reference; the dotted verticals are the structure's dated events. "
            "This rig is the ONLY window containing the build phase — both "
            "hindcast windows begin 15 years after installation, which is why "
            "neither of them can constrain f. The module reproduces the SIGN "
            "and the TIMING of the build and undershoots its AMPLITUDE: no "
            "admissible M matches the surveyed fillet on this grid, because a "
            "real fillet of about 190 m is being represented in a 500 m domain "
            "by a volume-neutral dipole whose real counterpart has a downdrift "
            "extent of 0 m. That gap is what the source/sink calibration "
            "absorbs. The observed peak is {opk:.0f} m in {pk}, the modelled "
            "peak {mpk:.0f} m, and at {endy} the survey reads {obs:+.0f} m "
            "against the model's {mod:+.0f} m. (b) The rate the module applied "
            "over the same years: full trapping while the structure was sound, "
            "the ramp from the 1996 last repair to the 2003 storm, then held "
            "at M·f. The rig runs 1967 off 1984 topography, deliberately."
            .format(M=applied_M, f=applied_f, opk=obs_vals.max(), pk=peak_year,
                    mpk=modelled.max(), endy=int(obs_years[-1]),
                    obs=final_obs, mod=final_mod))

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    written = save(figure, FIGURE_DIR / "full_life_1967_2017.png", close=True)
    print(f"wrote {written[0]}")


if __name__ == "__main__":
    main()
