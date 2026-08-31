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

Writes output/groin_sweep/figures/full_life_1967_2017.png

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import re
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

from HAT_groin_sweep_config import (  # noqa: E402
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

RAW_RUNS = PROJECT_BASE_DIR / "output" / "raw_runs"
GROIN_RUN = "HAT_1967_2018_edge_calibrated_groin"
NO_GROIN_RUN = "HAT_1967_2018_edge_calibrated_no_groin"


def _matrix(run_name):
    return RAW_RUNS / run_name / f"{run_name}_shoreline_matrix.npy"


SHORELINE = _matrix(GROIN_RUN)
NO_GROIN_SHORELINE = _matrix(NO_GROIN_RUN)
DIAGNOSTICS = (RAW_RUNS / GROIN_RUN
               / f"{GROIN_RUN}_groin_diagnostics.csv")

# A rig run directory does NOT name its own parameters -- the sweep writes every
# cell into one run name, so whatever survives is the last cell that finished.
# On 2026-08-30 this directory was found holding an UNSTABLE M = 70 cell whose
# fillet ran away to 444 m, while being named as though it were the calibrated
# run. So the applied rate is read back from the diagnostics and checked here
# rather than trusted from the label.
EXPECT_M, EXPECT_F = 60.0, 0.6

FIGURE_DIR = PROJECT_BASE_DIR / "output" / "groin_sweep" / "figures"

OBSERVED_COLOR, MODEL_COLOR, RATE_COLOR = "#1A1A1A", "#FF8C00", "#2A6FB0"
EVENT_COLOR, GRID_COLOR = "#C2185B", "#D5D8DD"

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

    plt.rcParams.update({"font.size": 10, "axes.linewidth": 0.7,
                         "legend.frameon": False, "pdf.fonttype": 42})
    figure, (ax_fillet, ax_rate) = plt.subplots(
        2, 1, figsize=(12.5, 8.2), sharex=True,
        gridspec_kw={"height_ratios": [2.4, 1.0], "hspace": 0.12})

    # ---- (a) fillet ------------------------------------------------------
    for year, label in EVENTS:
        ax_fillet.axvline(year, color=EVENT_COLOR, linewidth=1.0,
                          linestyle=(0, (4, 2)), alpha=0.75, zorder=2)
        ax_fillet.annotate(label, xy=(year, 1.005),
                           xycoords=("data", "axes fraction"),
                           rotation=90, ha="right", va="bottom",
                           fontsize=8.0, color=EVENT_COLOR)

    ax_fillet.axhline(0.0, color="#BBBBBB", linewidth=0.8, zorder=1)
    if no_groin is not None:
        ax_fillet.plot(years, no_groin, color="#777777", linestyle=":",
                       linewidth=1.9, label="model, groin OFF", zorder=3)
    ax_fillet.plot(years, modelled, color=MODEL_COLOR, linewidth=2.4,
                   label=f"model, rig at M = {applied_M:g}, f = {applied_f:.1f}",
                   zorder=4)
    ax_fillet.plot(obs_years, obs_vals, marker="o", markersize=6.5,
                   linestyle="none", color=OBSERVED_COLOR,
                   label=f"surveyed fillet ({len(obs_years)} dates)", zorder=5)

    peak_year = int(obs_years[int(np.argmax(obs_vals))])
    ax_fillet.set_ylabel("Fillet, D5 − D6 offset since 1967 (m)\n"
                         "[+ = updrift side holding]")
    ax_fillet.set_title(
        "The Buxton groin over its whole life, 1967–2017 — "
        "modelled fillet against every survey",
        fontsize=12.5, pad=26)
    ax_fillet.grid(axis="y", color=GRID_COLOR, linewidth=0.6)
    ax_fillet.set_axisbelow(True)
    ax_fillet.legend(loc="upper left", fontsize=9)
    for side in ("top", "right"):
        ax_fillet.spines[side].set_visible(False)

    final_obs = float(obs_vals[-1])
    final_mod = float(modelled[-1])
    ax_fillet.annotate(
        f"observed peak {obs_vals.max():.0f} m in {peak_year}\n"
        f"model peaks {modelled.max():.0f} m\n"
        f"end {int(obs_years[-1])}: observed {final_obs:+.0f} m, "
        f"model {final_mod:+.0f} m",
        xy=(0.985, 0.13), xycoords="axes fraction", ha="right", va="bottom",
        fontsize=8.8, color="#444444")

    # ---- (b) what the module applied -------------------------------------
    ax_rate.plot(diagnostics["model_year"],
                 diagnostics["trapping_rate_m_yr_applied"]
                 if "trapping_rate_m_yr_applied" in diagnostics
                 else diagnostics["trapping_rate_applied_m_yr"],
                 color=RATE_COLOR, linewidth=2.2, drawstyle="steps-post")
    for year, _ in EVENTS:
        ax_rate.axvline(year, color=EVENT_COLOR, linewidth=1.0,
                        linestyle=(0, (4, 2)), alpha=0.75, zorder=2)
    ax_rate.set_ylabel("Applied trapping\nrate M_eff (m/yr)")
    ax_rate.set_xlabel("Year")
    ax_rate.grid(axis="y", color=GRID_COLOR, linewidth=0.6)
    ax_rate.set_axisbelow(True)
    for side in ("top", "right"):
        ax_rate.spines[side].set_visible(False)
    ax_rate.annotate("M = 60 while sound        ramp 1996→2003        "
                     "hold at M·f = 36",
                     xy=(0.5, 0.06), xycoords="axes fraction", ha="center",
                     fontsize=8.8, color="#444444")

    figure.text(
        0.012, 0.012,
        "The rig is the ONLY window containing the build phase -- both hindcast "
        "windows begin 15 years after installation, which is why they cannot "
        "constrain f. The module reproduces the SIGN and TIMING of the build "
        "and undershoots its AMPLITUDE: no admissible M matches the fillet on "
        "this grid (a ~190 m real fillet in a 500 m domain, and a "
        "volume-neutral dipole where the real structure's downdrift extent is "
        "0 m). That gap is what the source/sink calibration absorbs. The rig "
        "runs 1967 off 1984 topography, deliberately.",
        fontsize=7.6, color="#555555", wrap=True, va="bottom")

    figure.subplots_adjust(top=0.90, bottom=0.155, left=0.085, right=0.985)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    out = FIGURE_DIR / "full_life_1967_2017.png"
    figure.savefig(out, dpi=200)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
