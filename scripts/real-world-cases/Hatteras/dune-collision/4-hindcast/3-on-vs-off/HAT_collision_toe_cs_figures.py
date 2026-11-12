"""
The calibration position-change figures of 2026-11-04 on the dune toe and Cs, in one colour scheme: the base or the
module off in grey, the step figures as one indigo-teal-gold gradient, the collision
module otherwise in purple (a plum-to-deep-purple ramp as Cs grows; where the two toes are compared, the 1996-97
toes dashed and the all-survey toes solid), the road rule in deep green, Rodanthe 90% in burnt orange.

    python scripts/real-world-cases/Hatteras/dune-collision/4-hindcast/3-on-vs-off/HAT_collision_toe_cs_figures.py [name ...]

Runs (HAT_dune_collision.py, all 1996-2009, Larson alone, dune growth 0.02-0.04 with the module, no road rule and
Rodanthe 40% unless the step adds them): cal_step1_base (the base, before any change), cal_step2_growth (module off),
cal_step3_module_la[_cs50|_cs100] (one toe per domain from the 1996 and 1997 surveys), cal_step3_module_la_alltoes
[_cs10|_cs30|_cs50|_cs100|_cs110] (one toe per domain averaged over every 1996-2009 survey, 1998 left out), and
cal_alltoes_cs110_road / _final (+ the 150 m road rule, + Rodanthe 90%).

Writes, under the workspace's dune-collision/4-hindcast/ (D:), the figures named in FIGURES.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-04
"""
from __future__ import annotations

import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_collision_on_off_position as P  # noqa: E402

sys.path.insert(0, str(P._HAT_DIR / "dune-collision" / "4-hindcast" / "1-calibration-1996-2009"))

# --- CONFIG ------------------------------------------------------------------
GREY = "#8C8C8C"                                   # the base / module off
PURPLE = "#652381"                                 # the collision module, one line
PURPLE_LIGHT = "#B877B8"                           # the module on the 1996-97 toes (dashed), set against PURPLE
RAMP = ("#D9A6CF", "#B877B8", "#8F4DA0", "#652381", "#3F0F5E")   # Cs 1e-5, 1.9e-5, 3e-5, 5e-5, 1e-4
GREEN = "#1B6E4F"                                  # the road rule
ORANGE = "#E07B39"                                 # Rodanthe 90%
# The step figures: one gradient, deep indigo -> teal -> gold (viridis), a shade further along for each change added
STEPS = ("#3B1F6B", "#21918C", "#E3B219")          # + module, + road rule, + Rodanthe 90%
DASH = (0, (5, 2.5))
BASE = ("cal_step1_base", "Base (before any changes)", GREY, 2.4)
OFF = ("cal_step2_growth", "Module off", GREY, 2.4)
ALL = "all-survey toes"
CAP = ("Larson et al. (2004) alone, dune growth 0.02-0.04, no road rule, villages 40%; one toe per domain from the "
       "1996 and 1997 surveys, or averaged over every 1996-2009 survey (1998 left out).")
def cs(mantissa, exponent):
    """The impact coefficient as a journal writes it, C_s = 1.9 x 10^-5 (matplotlib mathtext)."""
    return rf"$C_s$ = {mantissa} × 10$^{{{exponent}}}$"


def plain(label):
    """A legend label without the mathtext, for the caption."""
    return label.replace("$C_s$", "Cs").replace(" × 10$^{", "e").replace("}$", "")


TALL = set()   # figures to draw taller (height_scale 1.6); none: every figure the same size
START_TOES, ALL_TOES = "toes from the 1996–97 surveys", "toes from all 1996–2009 surveys"
FIGURES = {
    "1-module-on-off/base_vs_module_position_1996_2009.png": (
        "Shoreline change with and without the collision module, 1996–2009",
        (BASE, ("cal_step3_module_la", "Collision module (with dune growth 0.02–0.04)", PURPLE, 2.2))),
    "_superseded/automorph-toes/alltoes_position_1996_2009.png": (
        "Sensitivity of shoreline change to the dune-toe definition, 1996–2009",
        (OFF, ("cal_step3_module_la_alltoes", f"Collision module, {ALL_TOES}", PURPLE, 2.2),
         ("cal_step3_module_la", f"Collision module, {START_TOES} (current)", PURPLE_LIGHT, 2.0, DASH))),
    "_superseded/automorph-toes/domaintoes_cs_position_1996_2009.png": (
        "Sensitivity of shoreline change to the impact coefficient, start-survey toes, 1996–2009",
        (OFF, ("cal_step3_module_la", cs(1.9, -5), RAMP[1], 2.0),
         ("cal_step3_module_la_cs50", cs(5.0, -5), RAMP[3], 2.0),
         ("cal_step3_module_la_cs100", cs(1.0, -4), RAMP[4], 2.2))),
    "_superseded/automorph-toes/alltoes_cs_position_1996_2009.png": (
        "Sensitivity of shoreline change to the impact coefficient, 1996–2009",
        (OFF, ("cal_step3_module_la_alltoes", cs(1.9, -5), RAMP[1], 2.0),
         ("cal_step3_module_la_alltoes_cs50", cs(5.0, -5), RAMP[3], 2.0),
         ("cal_step3_module_la_alltoes_cs100", cs(1.0, -4), RAMP[4], 2.2))),
    "_superseded/automorph-toes/alltoes_cs_sweep_position_1996_2009.png": (
        "Sensitivity of shoreline change to the impact coefficient, 1996–2009",
        (OFF,) + tuple((v, cs(m, e), c, 2.0) for v, (m, e), c in zip(
            ("cal_step3_module_la_alltoes_cs10", "cal_step3_module_la_alltoes", "cal_step3_module_la_alltoes_cs30",
             "cal_step3_module_la_alltoes_cs50", "cal_step3_module_la_alltoes_cs100"),
            ((1.0, -5), (1.9, -5), (3.0, -5), (5.0, -5), (1.0, -4)), RAMP))),
    "2-model-steps/base_vs_cs110_position_1996_2009.png": (
        "Shoreline change with the collision module at the published impact coefficient, 1996–2009",
        (BASE, ("cal_step3_module_la_alltoes_cs110", f"Collision module, {cs(1.1, -4)}", PURPLE, 2.2))),
    "2-model-steps/cs110_road_position_1996_2009.png": (
        "Cumulative effect of the collision module and road rule, 1996–2009",
        (BASE, ("cal_step3_module_la_alltoes_cs110", f"+ collision module, {cs(1.1, -4)}", STEPS[0], 2.0),
         ("cal_alltoes_cs110_road", "+ 150 m road rule", STEPS[1], 2.2))),
    "2-model-steps/cs110_steps_position_1996_2009.png": (
        "Cumulative effect of the model changes, 1996–2009",
        (BASE, ("cal_step3_module_la_alltoes_cs110", f"+ collision module, {cs(1.1, -4)}", STEPS[0], 2.0),
         ("cal_alltoes_cs110_road", "+ 150 m road rule", STEPS[1], 2.0),
         ("cal_alltoes_cs110_final", "+ Rodanthe 90% overwash to the beach", STEPS[2], 2.2))),
}
# The step figure built up one line at a time for slides (same axes): level 1 the base alone ... level 4 every step
_STEPS = FIGURES["2-model-steps/cs110_steps_position_1996_2009.png"]
LEGEND_ALL = {}   # figure -> every run of its build-up, so each level keeps the full legend's space
for _k in range(1, len(_STEPS[1]) + 1):
    FIGURES[f"2-model-steps/slides/cs110_steps_level{_k}_position_1996_2009.png"] = (_STEPS[0], _STEPS[1][:_k])
    LEGEND_ALL[f"2-model-steps/slides/cs110_steps_level{_k}_position_1996_2009.png"] = _STEPS[1]
# -----------------------------------------------------------------------------


def main(names=None):
    for name, (title, runs) in FIGURES.items():
        if names and name not in names and Path(name).name not in names:
            continue
        lines = "; ".join(plain(r[1]) for r in runs)
        out, *_ = P.figure(1996, 2009, "HAT_dune_collision", runs, title, name,
                           f"each line one run ({lines}). {CAP}", height_scale=1.6 if name in TALL else 1.0,
                           legend_runs=LEGEND_ALL.get(name))
        print(f"wrote {out}")


if __name__ == "__main__":
    main(sys.argv[1:] or None)
