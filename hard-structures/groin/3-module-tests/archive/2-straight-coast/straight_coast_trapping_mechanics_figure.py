#!/usr/bin/env python3
"""
How the trapping groins work inside BRIE: the drift the waves carry, the exchange BRIE's diffusion
passes across a cell face, and what each trapping module blocks.

    python straight_coast_trapping_mechanics_figure.py   ->  comparison/figures/straight_coast_trapping_mechanics.png

(a) the wave climate's net drift Q(theta) and BRIE's diffusive exchange D(theta) tan(theta) depth on a
straight coast, against orientation; (b) the two, and what the trapping modules block, at three
orientations; (c) the bookkeeping at the groin on a bent coast, from the pinned run's tables;
(d) the 50-year sand totals of the pinned and conserving modules. Run straight_coast_test.py first.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-08
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "2-straight-coast" / "1-run"))   # the test, kept in place
import straight_coast_test as sc  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C, INK, INK_MUTED, _title, apply_style, figsize, open_frame, record_caption, save)

# --- CONFIG ------------------------------------------------------------------
OUT = Path(__file__).resolve().parent / "comparison" / "figures" / "straight_coast_trapping_mechanics.png"
EXAMPLES = (-30, 0, 20)          # orientations compared in (b)
BENT = (-30, 10)                 # orientation and year of the bent-coast bookkeeping in (c)
DRIFT_C, BRIE_C = INK, C["BASE"]
PIN_C, CONS_C, DRIFT_MOD_C = (sc.MODULES[k]["col"] for k in ("pinned", "conserving", "drift"))
# -----------------------------------------------------------------------------


def table(name, kind="by_year"):
    f = sc.MODULES[name]["folder"]
    return pd.read_csv(sc.module_dir(sc.MODULES[name]) / "tables" / f"straight_coast_{f}_{kind}.csv")


def bars(ax, groups, series, ylab):
    """Grouped bars: `series` is [(label, colour, values per group)]."""
    n = len(series)
    w = 0.8 / n
    x = np.arange(len(groups))
    for k, (lab, col, vals) in enumerate(series):
        ax.bar(x + (k - (n - 1) / 2) * w, vals, w * 0.92, color=col, label=lab, lw=0)
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.set_xticks(x, groups)
    ax.tick_params(axis="x", length=0)
    ax.set_ylabel(ylab)
    ax.grid(axis="y")
    open_frame(ax)


def main():
    apply_style()
    t = sc.brie_tables()
    th = np.arange(-40, 41)
    q0 = float(sc.drift(t, 0)[0])
    sgn = np.sign(q0)                                   # toward the groin is the drift's own direction at 0
    q = sgn * sc.drift(t, th)
    d = t["coast_diff"][sc.angle_index(th, t["climl"])]
    ex = sgn * d * np.tan(np.deg2rad(th)) * t["depth"]  # BRIE's diffusive exchange across a face, same sign
    k = np.flatnonzero(np.diff(np.sign(q)) != 0)
    zero_q = float(np.interp(0, q[k[0]:k[0] + 2][::-1], th[k[0]:k[0] + 2][::-1]))

    fig, axes = plt.subplots(2, 2, figsize=figsize("double", height=6.2), constrained_layout=True)

    # (a) what the waves carry vs what BRIE's diffusion passes
    ax = axes[0, 0]
    ax.axvspan(0, 40, color=C["BASE_FILL"], alpha=0.5, lw=0)
    ax.plot(th, q / 1e6, color=DRIFT_C, lw=1.8, label="net drift the waves carry, Q($\\theta$)")
    ax.plot(th, ex / 1e6, color=BRIE_C, lw=1.6, ls="--", label="what BRIE's diffusion passes, D($\\theta$) tan$\\theta$ h")
    ax.axhline(abs(q0) / 1e6, color=DRIFT_C, lw=0.6, ls=":")
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.text(-39, abs(q0) / 1e6 - 0.06, f"Q at 0$^\\circ$: {abs(q0) / 1e6:.2f} M m$^3$/yr", fontsize=6.5, color=INK_MUTED, va="top")
    ax.text(20, 1.52, "BRIE's exchange runs\nagainst the drift", ha="center", fontsize=6.5, color=INK_MUTED, va="top")
    for e in EXAMPLES:
        ax.axvline(e, color=INK_MUTED, lw=0.4, ls=(0, (1, 2)))
    ax.set_ylim(-2.0, 1.6)
    ax.set_xlim(-40, 40)
    ax.set_xlabel("Shoreline orientation ($^\\circ$)")
    ax.set_ylabel("Sand carried toward the groin\n(M m$^3$/yr; negative = away)")
    ax.legend(loc="lower left", fontsize=6.5, frameon=False)
    ax.grid(axis="y")
    open_frame(ax)
    _title(ax, 0, "What moves the sand?")

    # (b) at three orientations: the two flows, and what each trapping module blocks
    ax = axes[0, 1]
    qe = np.interp(EXAMPLES, th, q) / 1e6
    ee = np.interp(EXAMPLES, th, ex) / 1e6
    bars(ax, [f"{e:+d}$^\\circ$" for e in EXAMPLES], [
        ("net drift Q", DRIFT_C, qe),
        ("BRIE's exchange", BRIE_C, ee),
        (f"pinned / conserving block {sc.BLOCKING_FRACTION} $\\times$ exchange", PIN_C, sc.BLOCKING_FRACTION * ee),
        (f"drift version blocks {sc.BLOCKING_FRACTION} $\\times$ Q", DRIFT_MOD_C, sc.BLOCKING_FRACTION * qe),
    ], "Toward the groin (M m$^3$/yr)")
    ax.set_xlabel("Shoreline orientation")
    ax.set_ylim(-1.05, 1.15)
    ax.legend(loc="lower left", fontsize=6.3, frameon=False)
    _title(ax, 1, "What does each module block?")

    # (c) the bookkeeping at the groin once the coast bends, from the pinned run
    th_b, yr_b = BENT
    pin = table("pinned")
    row = pin[(pin.theta0_deg == th_b) & (pin.year == yr_b)].iloc[0]
    s = abs(np.tan(np.deg2rad(row.face_angle_deg))) * sc.DY_M          # the step across the groin, m
    f0 = float(sc.drift(t, th_b)[0])
    r_up, r_down = (row.r_hi, row.r_lo) if f0 < 0 else (row.r_lo, row.r_hi)
    r_face = row.r_lo
    brie = np.array([2 * r_down * s, 2 * r_up * s])
    ax = axes[1, 0]
    bars(ax, ["downdrift cell", "updrift cell"], [
        ("BRIE's exchange, each cell's own r", BRIE_C, brie),
        ("pinned cancels (own r)", PIN_C, sc.BLOCKING_FRACTION * brie),
        ("conserving cancels (the face's r)", CONS_C, sc.BLOCKING_FRACTION * np.array([2 * r_face * s] * 2)),
    ], "Shoreline moved at the groin\n(m/yr)")
    ax.set_ylim(0, brie.max() * 1.45)
    for xi, v in enumerate(sc.BLOCKING_FRACTION * brie):
        ax.text(xi - 0.0, v + brie.max() * 0.03, f"{v:.0f}", ha="center", fontsize=6.5, color=PIN_C)
    ax.text(0.98, 0.72, f"coast at {th_b:+d}$^\\circ$, year {yr_b}\nstep {s:.0f} m, r$_{{down}}$ {r_down:.3f}, "
            f"r$_{{up}}$ {r_up:.3f}", transform=ax.transAxes, ha="right", va="top", fontsize=6.3, color=INK_MUTED)
    ax.legend(loc="upper left", fontsize=6.3, frameon=False)
    _title(ax, 2, "When the coast bends")

    # (d) where the sand goes in 50 years: module, solve, total
    summ = pd.read_csv(sc.HERE / "tables" / "straight_coast_summary_year50.csv")
    sm = summ[summ.theta0_deg == th_b].set_index("module")
    ax = axes[1, 1]
    names = ["pinned", "conserving"]
    bars(ax, [sc.MODULES[n]["label"] for n in names], [
        ("the module adds", INK_MUTED, [sm.loc[n, "applied_net_cum_m"] for n in names]),
        ("BRIE's solve adds", BRIE_C, [sm.loc[n, "solve_error_m"] for n in names]),
        ("total", INK, [sm.loc[n, "reach_change_m"] for n in names]),
    ], "Sand gained in 50 yr\n(m of summed shoreline)")
    for xi, n in enumerate(names):
        ax.text(xi + 0.27, sm.loc[n, "reach_change_m"] + 60, f"{sm.loc[n, 'reach_change_m']:+.0f}",
                ha="center", va="bottom", fontsize=6.5, color=INK)
    ax.text(0.98, 0.03, f"coast at {th_b:+d}$^\\circ$", transform=ax.transAxes, ha="right", fontsize=6.3,
            color=INK_MUTED)
    ax.legend(loc="upper left", fontsize=6.3, frameon=False)
    _title(ax, 3, "Which one conserves sand?")
    fig.suptitle("How do the trapping groins work inside BRIE?")
    save(fig, OUT, dpi=250, close=True)

    i30, i20 = (EXAMPLES.index(e) for e in (-30, 20))
    record_caption(OUT, (
        "**Tests:** nothing new; this figure explains the mechanics behind the trapping groin modules, "
        "why the pinned and conserving versions trap on the wrong side on some coasts, and why only the "
        "pinned version conserves sand. **How to read:** (a) On a straight coast at each orientation, the "
        "net longshore drift the option A wave climate carries toward the groin (solid; BRIE's own "
        "_coast_qs convolved with the wave-angle pdf) and the exchange BRIE's diffusion passes across a cell "
        f"face, D($\\theta$) tan($\\theta$) times the {t['depth']:.1f} m active depth (dashed). BRIE moves the "
        "shoreline by diffusion only, with D = -(1/depth) dQ/d$\\theta$, so a drift that is the same in "
        f"every cell, the {abs(q0) / 1e6:.2f} M m$^3$/yr at 0 degrees (dotted), moves nothing and never "
        "appears; the exchange is zero at 0 degrees, follows the drift's direction on coasts turned against "
        f"it and runs against it on coasts turned with it (shaded). The drift reaches zero at {zero_q:+.0f} "
        "degrees. (b) At three orientations: the drift, BRIE's exchange, what the pinned and conserving "
        f"modules block ({sc.BLOCKING_FRACTION} of the exchange) and what the drift module blocks ({sc.BLOCKING_FRACTION} of the drift); "
        "positive is toward the groin, so a positive block builds the updrift side. (c) The bookkeeping at "
        f"the groin on a bent coast: the pinned run at {th_b:+d} degrees in year {yr_b}, when the fillet has "
        "steepened the step across the groin. BRIE's solve is row-scaled: each cell's row uses that cell's "
        "own diffusivity r for both its neighbours, so the exchange is counted once with the downdrift r and "
        "once with the updrift r, and the two counts differ. The pinned module cancels a share of each count "
        "as BRIE makes it; the conserving module uses the groin face's r for both sides. (d) The 50-year "
        f"sand totals of the two at {th_b:+d} degrees: what the module adds, what BRIE's solve adds, and "
        "the sum. **Shows:** the pinned and conserving modules block a share of BRIE's diffusive exchange, "
        "not of the drift: "
        f"{sc.BLOCKING_FRACTION * ee[i30]:+.2f} M m$^3$/yr at -30 degrees where the drift is {qe[i30]:.2f}, nothing at 0 "
        f"degrees where it is {qe[1]:.2f}, and {sc.BLOCKING_FRACTION * ee[i20]:+.2f} at +20 degrees, where blocking an "
        "exchange that runs against the drift builds the downdrift side. The drift module blocks "
        f"{sc.BLOCKING_FRACTION * qe[i30]:.2f}, {sc.BLOCKING_FRACTION * qe[1]:.2f} and {sc.BLOCKING_FRACTION * qe[i20]:.2f}. Once the coast bends the two "
        f"sides' r differ ({r_down:.3f} vs {r_up:.3f}), so the pinned module takes "
        f"{sc.BLOCKING_FRACTION * brie[0]:.0f} m from the downdrift cell and gives {sc.BLOCKING_FRACTION * brie[1]:.0f} m to the updrift "
        "cell; that mismatch mirrors BRIE's own, and over 50 years the solve hands it back "
        f"({sm.loc['pinned', 'applied_net_cum_m']:+.0f} m and {sm.loc['pinned', 'solve_error_m']:+.0f} m, total "
        f"{sm.loc['pinned', 'reach_change_m']:+.0f} m). The conserving module balances its own books, but "
        "BRIE's mismatch is then left in place and the reach gains "
        f"{sm.loc['conserving', 'reach_change_m']:+.0f} m."))
    print(f"wrote {OUT.relative_to(sc.REPO)}")


if __name__ == "__main__":
    main()
