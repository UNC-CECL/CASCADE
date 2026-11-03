"""
The dune collision module over 200 years: do the dunes level off where growth and collision balance, or keep
falling to the toe?

    python collision_long_term_test.py [--years 200]
        ->  figures/collision_long_term.png, figures/supporting/long_term_summary.csv

Barrier3D's own default test domain (tests/test_coupling, prefix barrier3d-default: 50 cells, berm 1.44 m and
dune ceiling 2.94 m MHW, its 1000-year Virginia Coast Reserve storm series), run alone for --years with the
module off, on with no beach, and on behind a 25 m beach (W0 45 m), at Cs 1.24e-4: A, the domain as it is
(own growth 0.55-0.95, toe 0.1 m above the berm); B, Hatteras-like (ceiling 6 m and dunes starting 4 m above
the berm, growth 0.02-0.04, toe 0.6 m above the berm, a Hatteras toe); C, B with A's low toe. Each
year is seeded as tests/test_dune_collision.py seeds it, so the runs share their random draws. First the
module runs on the series as it is, to see whether the module's stated limits stop it; the runs then cap the
storm periods at the module's 20 s. Dune migration with the shoreline is switched off in every run, so the dunes
stay put and only growth, overwash and collision act on them. Takes about two minutes.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-03
"""
from __future__ import annotations

import argparse
import shutil
import sys
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import yaml  # noqa: E402

_HERE = Path(__file__).resolve()
_HAT_DIR = next(_p for _p in _HERE.parents if (_p / "hatteras_paths.py").exists())
_REPO = next(_p for _p in _HERE.parents if (_p / "pyproject.toml").exists())
for _p in (_HAT_DIR, _REPO):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))
from hatteras_paths import workspace_twin  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C, INK_MUTED, _title, apply_style, caption, figsize, open_frame, save, support_dir)
from barrier3d import Barrier3d  # noqa: E402
from cascade.dune_collision import DuneCollisionManager, DuneCollisionOutsideTestedRange  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
DATA_DIR = _REPO / "tests" / "test_coupling"
PREFIX = "barrier3d-default"
YEARS = 200
CS = 1.24e-4
TOE_ABOVE_BERM_M = {"A": 0.1, "B": 0.6}  # m; B at a Hatteras toe (2.4 m NAVD88 = 2.0 m MHW)
BEACH_M = 25.0                   # the beach for the "behind a beach" run
W0_M = 45.0
# The two domains: the test domain as it is (low dunes, fast growth: overwash rules), and Hatteras-like (tall
# dunes, slow growth: the storms mostly reach the dunes in collision). Overrides of the test domain's inputs, m
GROWTH = {"A: test domain (low dunes, fast growth)": {"toe_above_berm_m": TOE_ABOVE_BERM_M["A"]},
          "B: Hatteras-like (tall dunes, slow growth)": {"Dmaxel": 1.9 + 6.0, "start_dune_m": 4.0,
                                                         "rmin": 0.02, "rmax": 0.04,
                                                         "toe_above_berm_m": TOE_ABOVE_BERM_M["B"]},
          "C: as B, toe 0.1 m above the berm": {"Dmaxel": 1.9 + 6.0, "start_dune_m": 4.0, "rmin": 0.02,
                                                "rmax": 0.04, "toe_above_berm_m": TOE_ABOVE_BERM_M["A"]}}
RUNS = {"module off": None, "module on, no beach": None, f"module on, {BEACH_M:.0f} m beach": BEACH_M}
SETTLE_WINDOW = (100, 190)      # years the equilibrium is the mean over: before the series' largest storm
BIG_STORM_YEAR = 193            # the 1000-year series' largest storm (4.04 m MHW for 120 h)
SETTLED_WITHIN_M = 0.1          # a run has settled when its 20-year mean first comes this close to that
OUT = workspace_twin(_HERE.parent) / "figures"
COLORS = {"module off": C["BASE"], "module on, no beach": C["ACCENT"],
          f"module on, {BEACH_M:.0f} m beach": C["ACCENT_FILL"]}
# -----------------------------------------------------------------------------


# The test domain's inputs, copied with a longer run and, optionally, other growth bounds
def domain_inputs(years, overrides):
    d = Path(tempfile.mkdtemp(prefix="collision_long_term_"))
    for f in DATA_DIR.iterdir():
        if f.is_file():
            shutil.copy(f, d)
    yml = d / f"{PREFIX}-parameters.yaml"
    p = yaml.safe_load(yml.read_text())
    p["TMAX"] = years + 1
    p.update({k: v for k, v in overrides.items() if k not in ("start_dune_m", "toe_above_berm_m")})
    yml.write_text(yaml.safe_dump(p))
    return d


# One run: Barrier3D alone, the module as given; mean crest above the berm and collision lowering each year
def run(years, overrides, beach_m, module, cap_period=True):
    folder = domain_inputs(years, overrides)
    np.random.seed(0)
    b3d = Barrier3d.from_path(folder, prefix=PREFIX)
    b3d._dune_migration_on = False     # the dune line stays put: only growth, overwash and collision act
    if "start_dune_m" in overrides:    # the test domain reads its start dunes from a file: set them here
        b3d.DuneDomain[0, :, :] = overrides["start_dune_m"] / 10
    if cap_period:
        b3d._StormSeries[:, 3] = np.minimum(b3d._StormSeries[:, 3], DuneCollisionManager.STORM_PERIOD_RANGE_S[1])
    n = b3d.BarrierLength
    mgr = None
    if module:
        toe = np.full(n, b3d._BermEl + overrides["toe_above_berm_m"] / 10)
        width = None if beach_m is None else np.full(n, beach_m / 10)
        mgr = DuneCollisionManager(b3d, toe, cs=CS, beach_width=width, beach_width_limit=W0_M / 10)
    stopped = None
    try:
        for _ in range(years):
            np.random.seed(1000 + b3d.time_index)
            b3d.update()
            b3d.update_dune_domain()
            if b3d.drown_break:
                stopped = f"drowned at year {b3d.time_index}"
                break
    except DuneCollisionOutsideTestedRange as e:
        stopped = f"year {b3d.time_index}: {e}"
    shutil.rmtree(folder, ignore_errors=True)
    t = b3d.time_index
    crest = b3d.DuneDomain[:t].max(axis=2).mean(axis=1) * 10               # m above the berm
    low = np.zeros(t) if mgr is None else np.r_[0, np.array(mgr.crest_lowering_TS) * 10][:t]
    return dict(crest=crest, low=np.pad(low, (0, t - len(low))), stopped=stopped,
                dmax=(b3d._Dmaxel - b3d._BermEl) * 10)


# The equilibrium (mean over SETTLE_WINDOW), its drift (the window's second half less its first), and the
# first year the 20-year mean comes within SETTLED_WITHIN_M of it
def settle(crest):
    lo, hi = SETTLE_WINDOW
    mid = (lo + hi) // 2
    eq = float(np.mean(crest[lo:hi]))
    drift = float(np.mean(crest[mid:hi]) - np.mean(crest[lo:mid]))
    roll = pd.Series(crest).rolling(20).mean().to_numpy()
    year = next((i for i in range(len(roll)) if np.abs(roll[i] - eq) <= SETTLED_WITHIN_M), None)
    return eq, drift, year


# Run: the unclipped check, the six runs, the table, the figure
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--years", type=int, default=YEARS)
    args = ap.parse_args()

    # Does the module's tested range stop it on the series as it is?
    raw = {g: run(args.years, o, None, True, cap_period=False)["stopped"] or "ran to the end"
           for g, o in GROWTH.items()}
    print(f"series as it is, module on: {raw}")

    rows, res = [], {}
    for gname, growth in GROWTH.items():
        for rname, beach in RUNS.items():
            r = run(args.years, growth, beach, rname != "module off")
            eq, drift, yr = settle(r["crest"])
            res[gname, rname] = r
            rows.append(dict(growth=gname, run=rname, equilibrium_crest_m_above_berm=round(eq, 3),
                             drift_window_halves_m=round(drift, 3), settled_by_year=yr, dune_ceiling_m_above_berm=round(r["dmax"], 2),
                             mean_collision_lowering_m_per_yr=round(float(np.mean(r["low"])), 4),
                             crest_after_big_storm_m=round(float(r["crest"][min(BIG_STORM_YEAR + 1, len(r["crest"]) - 1)]), 3),
                             crest_at_end_m=round(float(r["crest"][-1]), 3),
                             share_years_crest_below_0p2m=round(float(np.mean(r["crest"] < 0.2)), 3),
                             stopped=r["stopped"] or ""))
            print(f"{gname:34s} {rname:24s} equilibrium {eq:.2f} m above berm, drift {drift:+.2f} m, within 0.1 m by year {yr}, "
                  f"collision {np.mean(r['low']):.3f} m/yr" + (f"  [{r['stopped']}]" if r["stopped"] else ""))
    tab = pd.DataFrame(rows)
    sup = support_dir(OUT)
    tab.to_csv(sup / "long_term_summary.csv", index=False)

    apply_style()
    fig, axes = plt.subplots(len(GROWTH), 2, figsize=figsize("double", 0.86), layout="constrained", sharex=True,
                             gridspec_kw=dict(width_ratios=(1.6, 1)))
    for i, gname in enumerate(GROWTH):
        a, b = axes[i]
        for rname in RUNS:
            r = res[gname, rname]
            yrs = np.arange(len(r["crest"]))
            a.plot(yrs, pd.Series(r["crest"]).rolling(5, min_periods=1).mean(), color=COLORS[rname], lw=1.3,
                   label=rname)
            if rname != "module off":
                b.plot(yrs, pd.Series(r["low"]).rolling(5, min_periods=1).mean(), color=COLORS[rname], lw=1.3)
        dmax = res[gname, "module off"]["dmax"]
        a.axhline(dmax, color=INK_MUTED, lw=0.7, ls=":")
        a.axvspan(*SETTLE_WINDOW, color="0.95", lw=0, zorder=0)
        a.axvline(BIG_STORM_YEAR, color=C["LATE"], lw=0.7, ls=":", zorder=0)
        toe_m = GROWTH[gname]["toe_above_berm_m"]
        a.axhline(toe_m, color=INK_MUTED, lw=0.7, ls="--")
        a.text(args.years, dmax + 0.03, "dune ceiling", ha="right", va="bottom", fontsize=7, color=INK_MUTED)
        a.text(args.years, toe_m + 0.03, "toe", ha="right", va="bottom", fontsize=7, color=INK_MUTED)
        a.set(ylim=(0, dmax * 1.15), ylabel="Mean crest above\nthe berm (m)")
        b.set(ylabel="Collision lowering\n(m/yr)")
        _title(a, 2 * i, gname)
        _title(b, 2 * i + 1, "Collision lowering")
        for ax in (a, b):
            open_frame(ax)
    for ax in axes[-1]:
        ax.set_xlabel("Model year")
    h, lab = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, lab, loc="outside upper center", ncol=3, frameon=False, fontsize=8)

    def fmt(g, r):
        x = tab[(tab.growth == g) & (tab.run == r)].iloc[0]
        return (f"{x.equilibrium_crest_m_above_berm:.2f} m (drift {x.drift_window_halves_m:+.2f} m between the window's halves, "
                f"within 0.1 m by year {x.settled_by_year})")
    g1, g2, g3 = list(GROWTH)
    on_beach = f"module on, {BEACH_M:.0f} m beach"
    caption(fig, f"**Tests:** the dune collision module on its own over {args.years} years, no Hatteras fit: do the "
            f"dunes level off where growth and collision balance, or fall to the toe? **How to read:** Barrier3D's "
            f"default test domain (tests/test_coupling: 50 cells, its Virginia Coast Reserve storm series, peak water "
            f"up to 4.0 m above MHW), run alone: (a, b) as it is, dune ceiling "
            f"{res[g1, 'module off']['dmax']:.2f} m above the berm, dunes starting 0.5 m above it, its own growth "
            f"0.55-0.95; (c, d) Hatteras-like, ceiling {res[g2, 'module off']['dmax']:.2f} m, dunes starting 4 m "
            f"above the berm, growth 0.02-0.04. Each with the module off, on with no beach, and on behind a "
            f"{BEACH_M:.0f} m beach (W0 {W0_M:.0f} m), Cs {CS:g}, toe {TOE_ABOVE_BERM_M['A']} m (A) and {TOE_ABOVE_BERM_M['B']} m (B, a Hatteras toe; C, B with A's toe) above the "
            f"berm. Dune migration with the shoreline is off, so only growth, overwash and collision act. Left, "
            f"the mean crest above the berm (5-year running mean); right, the module's crest lowering each year. "
            f"Every year is seeded alike, so the runs share their storms. The storm periods are capped at the "
            f"module's 20 s (on the series as it is, module on: "
            + "; ".join(f"{g.split(':')[0]} {s}" for g, s in raw.items())
            + f"). **Shows:** the mean crest over years {SETTLE_WINDOW[0]}-{SETTLE_WINDOW[1]} (before the series' largest "
            f"storm, year {BIG_STORM_YEAR}, 4.04 m for 120 h, which every run feels at the end): on A, off {fmt(g1, 'module off')}, "
            f"on {fmt(g1, 'module on, no beach')}, behind the beach {fmt(g1, on_beach)}; on B, off "
            f"{fmt(g2, 'module off')}, on {fmt(g2, 'module on, no beach')}, behind the beach {fmt(g2, on_beach)}; on C "
            f"(B with the low toe), off {fmt(g3, 'module off')}, on {fmt(g3, 'module on, no beach')}, behind the "
            f"beach {fmt(g3, on_beach)}. The balance depends on the toe: Barrier3D's logistic dune growth is "
            f"slowest on low dunes, so once collision outpaces it the dunes fall to the toe and do not come back.")
    save(fig, OUT / "collision_long_term", close=True)
    print(tab.to_string(index=False))


if __name__ == "__main__":
    main()
