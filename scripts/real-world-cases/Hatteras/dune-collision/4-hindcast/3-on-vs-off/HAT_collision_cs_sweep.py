"""
The collision module's Cs, in the model: Larson alone (no beach factor, the default since 2026-11-04) at five Cs,
calibration (1996-2009) and test (2009-2025), against the base with the module off.

    python scripts/real-world-cases/Hatteras/dune-collision/4-hindcast/3-on-vs-off/HAT_collision_cs_sweep.py

Cs 1e-5, 1.9e-5 (the offline fit, larson_storm in the rule comparison), 3e-5, 5e-5 and 8e-5, each with one toe per
domain, the 198 m road rule and Rodanthe 90%, as the default hindcast. Cs is chosen on the calibration; the
test runs show how much the choice matters. Also Larson x the domain beach factor (W0 25 m, Cs 1.02e-4, the default
until 2026-11-04) for comparison.

Runs: HAT_dune_collision.py cal_la_cs10, cal_step6_toeavg_la (1.9e-5), cal_la_cs30, cal_la_cs50, cal_la_cs80,
cal_step6_toeavg (beach factor), cal_step1_base; HAT_collision_test_2009.py la_cs10, nocape_rd198_dry_toeavg_la,
la_cs30, la_cs50, la_cs80, nocape_rd198_dry_toeavg, step1_base.

Writes, in the workspace's dune-collision/4-hindcast/4-cs-and-toes/ (D:): cs_sweep_fit.png (shoreline and crest error
against Cs, both periods), cs_sweep_position_<start>_<end>.png (shoreline position change at each Cs, CoastSat)
and supporting/cs_sweep_fit.csv.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-04
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_collision_model_steps as M  # noqa: E402
import HAT_collision_on_off_position as P  # noqa: E402
OUT = P.OUT / "4-cs-and-toes"   # this script's figures (topic subfolder)

# --- CONFIG ------------------------------------------------------------------
CS_VALUES = (1e-5, 1.9e-5, 3e-5, 5e-5, 8e-5)
OFFLINE_CS = 1.9e-5
PERIODS = {
    (1996, 2009): ("HAT_dune_collision", "Calibration 1996–2009", "cal_step1_base", "cal_step6_toeavg",
                   ("cal_la_cs10", "cal_step6_toeavg_la", "cal_la_cs30", "cal_la_cs50", "cal_la_cs80")),
    (2009, 2025): ("HAT_collision_test_2009", "Test 2009–2025", "step1_base", "nocape_rd198_dry_toeavg",
                   ("la_cs10", "nocape_rd198_dry_toeavg_la", "la_cs30", "la_cs50", "la_cs80")),
}
# Cs on one sequential purple ramp, light (weak) to dark (strong); not blue, which is CoastSat's
CS_COLORS = ("#BCBDDC", "#9E9AC8", "#807DBA", "#6A51A3", "#3F007D")
BASE_C, BEACH_C = "#7F7F7F", "#D95F02"
# -----------------------------------------------------------------------------


def cs_label(cs):
    m, e = f"{cs:.1e}".split("e")
    return f"Cs {float(m):g}e{int(e)}"


def fit_table():
    rows = []
    for (start, end), (module, _, base, beach, sweep) in PERIODS.items():
        steps = [("base", base), ("beach", beach)] + [(f"cs{cs:g}", v) for cs, v in zip(CS_VALUES, sweep)]
        t = M.fit(module, steps, start)
        t["cs"] = t.step.map({f"cs{cs:g}": cs for cs in CS_VALUES})
        rows.append(t)
    return pd.concat(rows, ignore_index=True)


def fit_figure(t):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer.hat_figure_style import record_caption
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7), sharex=True, sharey="row", constrained_layout=True)
    for j, ((start, end), (_, title, *_)) in enumerate(PERIODS.items()):
        d = t[(t.start == start) & (t.zone == "island")]
        for i, (col, ylab) in enumerate((("shoreline_rmse_m_yr", "Shoreline change rate RMSE (m/yr)"),
                                         ("crest_rmse_m", "Dune crest RMSE against lidar (m)"))):
            ax = axes[i, j]
            s = d.dropna(subset=["cs"]).sort_values("cs")
            ax.plot(s.cs, s[col], "-o", color=CS_COLORS[3], lw=2, ms=6, label="Larson alone", zorder=4)
            for x, y in zip(s.cs, s[col]):
                ax.annotate(f"{y:.2f}", (x, y), textcoords="offset points", xytext=(0, 8), ha="center",
                            fontsize=8, color="#2B2B2B")
            for step, c, ls, lab in (("base", BASE_C, "--", "Module off (base)"),
                                     ("beach", BEACH_C, ":", "Larson × beach factor, Cs 1.02e-4")):
                y = d[d.step == step][col].iloc[0]
                ax.axhline(y, color=c, ls=ls, lw=1.6, label=f"{lab} ({y:.2f})" if i == 0 and j == 0 else lab)
            ax.axvline(OFFLINE_CS, color="#BDBDBD", lw=0.8, zorder=0)
            ax.set_xscale("log")
            ax.set_xticks(CS_VALUES, [f"{cs:.1e}".replace("e-0", "e-") for cs in CS_VALUES], fontsize=8.5)
            ax.minorticks_off()
            ax.grid(axis="y", color="#E0E0E0", lw=0.8)
            ax.spines[["top", "right"]].set_visible(False)
            if j == 0:
                ax.set_ylabel(ylab)
            if i == 0:
                ax.set_title(title, fontsize=12, fontweight="bold")
            if i == 1:
                ax.set_xlabel("Cs (log scale; grey line, the offline fit)")
    h, lab = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, ["Larson alone", "Module off (base)", "Larson × beach factor, Cs 1.02e-4"], loc="lower center",
               bbox_to_anchor=(0.5, -0.06), ncol=3, frameon=False, fontsize=9.5)
    fig.suptitle("Collision module, Larson alone: model fit against Cs (lower is better)", fontsize=13,
                 fontweight="bold")
    out = OUT / "cs_sweep_fit.png"
    fig.savefig(out, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    record_caption(out, (
        "Model fit against the collision module's Cs, Larson alone (no beach factor), island GIS 6-90; lower is "
        "better. Top: RMSE of the modelled shoreline change rate against CoastSat's; bottom: RMSE of the modelled "
        "domain-median dune crest against the lidar, every survey in the period. Every run with one toe per domain, "
        "the 198 m road rule and Rodanthe 90%. Dashed, the base (module off, no road rule, Rodanthe 40%); dotted, "
        "the module with the domain beach factor (W0 25 m, Cs 1.02e-4). Grey line: Cs 1.9e-5, the offline fit. Cs is "
        "chosen on the calibration (left); the test (right) shows how much the choice matters."))
    return out


def main():
    t = fit_table()
    (OUT / "supporting").mkdir(parents=True, exist_ok=True)
    t.round(4).to_csv(OUT / "supporting" / "cs_sweep_fit.csv", index=False)
    print(f"wrote {fit_figure(t)}")
    for (start, end), (module, title, base, _, sweep) in PERIODS.items():
        runs = ((base, "Base (module off)", BASE_C, 2.2),) + tuple(
            (v, f"Larson alone, {cs_label(cs)}" + (" (offline fit)" if cs == OFFLINE_CS else ""), c,
             2.4 if cs == OFFLINE_CS else 1.6)
            for cs, v, c in zip(CS_VALUES, sweep, CS_COLORS))
        out, *_ = P.figure(
            start, end, module, runs, f"Collision module, Larson alone, at five Cs, {start}–{end}",
            f"4-cs-and-toes/cs_sweep_position_{start}_{end}.png",
            "the base (module off, no 198 m road rule, Rodanthe 40%) and the collision module with Larson alone (no "
            "beach factor) at Cs 1e-5 to 8e-5, light to dark, each with one toe per domain, the 198 m road rule and "
            "Rodanthe 90%; Cs 1.9e-5 is the offline fit.")
        print(f"wrote {out}")
        # The two candidates alone: the in-model best (1e-5) and the offline fit (1.9e-5)
        pair = ((base, "Base (module off)", BASE_C, 2.2),
                (sweep[0], "Larson alone, Cs 1e-5 (in-model best, calibration)", "#E7298A", 2.0),
                (sweep[1], "Larson alone, Cs 1.9e-5 (offline fit)", "#6A51A3", 2.2))
        out, *_ = P.figure(
            start, end, module, pair, f"Collision module, Larson alone: Cs 1e-5 against 1.9e-5, {start}–{end}",
            f"4-cs-and-toes/cs_pair_position_{start}_{end}.png",
            "the base (module off, no 198 m road rule, Rodanthe 40%) and the collision module with Larson alone (no "
            "beach factor) at Cs 1e-5 (the best calibration fit in the model) and 1.9e-5 (the offline fit to the lidar "
            "crest changes), each with one toe per domain, the 198 m road rule and Rodanthe 90%.")
        print(f"wrote {out}")
    pd.set_option("display.width", 200)
    for col in ("shoreline_rmse_m_yr", "crest_rmse_m", "crest_bias_m"):
        print(f"\n{col}")
        print(t.pivot_table(index=["start", "zone"], columns="step", values=col, sort=False).round(3).to_string())


if __name__ == "__main__":
    main()
