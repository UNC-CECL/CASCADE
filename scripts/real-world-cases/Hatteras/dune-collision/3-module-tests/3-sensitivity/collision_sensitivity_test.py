"""
How much the dune collision module's crest loss depends on its two parameters, Cs and the beach limit W0.

    python collision_sensitivity_test.py [--gis 12] [--year 1999]
        ->  figures/collision_sensitivity.png, figures/supporting/sensitivity.csv

Offline, no model runs: the domain's 50 dune cells start the year at their lidar crest (the survey before
the year), with the toes and beach widths the hindcast reads, and each storm of the year in the hindcast
storm series lowers them through the module's own functions, exactly as
../../2-module-explained/collision_year_gif.py does (its cells() and run_year(), imported unchanged; only
the Cs and W0 they read are set here). Cs from 0 to the module's tested maximum, 2.5e-4, at W0 of 10-60 m (dry beach, start to toe);
plus Larson et al.'s (2004) laboratory value, 1.4e-3, outside the module's bounds (the functions take it;
DuneCollisionManager would stop a run with it). The lidar's own crest change over the year is the reference.
Takes seconds.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-03
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parents[2] / "2-module-explained"))
import collision_explained_figures as E  # noqa: E402  (paths, units, the storm series; Cs and W0 set below)
import collision_year_gif as G  # noqa: E402  (the cells and the year's storms through the module)
from cascade.dune_collision import DuneCollisionManager  # noqa: E402
from hatteras_paths import workspace_twin  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C, INK_MUTED, _title, apply_style, caption, figsize, open_frame, save, support_dir)

# --- CONFIG ------------------------------------------------------------------
GIS = 12
YEAR = 1999
CS_VALUES = np.linspace(0, DuneCollisionManager.CS_MAX, 26)
W0_VALUES = (10.0, 15.0, 25.0, 40.0, 60.0)     # m, dry beach (start to toe)
CS_CALIBRATED, W0_CALIBRATED = 1.02e-4, 25.0     # hat_run.yaml since 2026-11-03
CS_LARSON_FIELD, CS_LARSON_LAB = 1.7e-4, 1.4e-3  # Larson et al. (2004): US East Coast storms; laboratory
CS_LINES = (0.5e-4, CS_CALIBRATED, CS_LARSON_FIELD, DuneCollisionManager.CS_MAX)
OUT = workspace_twin(_HERE.parent) / "figures"
W0_RAMP = ("#d9c4e3", "#c2a5cf", "#9970ab", "#7b3294", "#40004b")
# -----------------------------------------------------------------------------


# The year's mean crest lowering over the domain's cells (m) at one Cs and W0
def lowering(d, st, cs, w0):
    E.CS, E.W0 = cs, w0
    steps = G.run_year(d, st)
    return abs(float(np.nanmean(d.crest_m_navd88.to_numpy() - steps[-1]["after"]))) + 0.0


# The lidar's own mean crest change over the year, on the same cells (the survey after less the one before)
def lidar_change(d, year, gis):
    c = pd.read_csv(G.CRESTS)
    c = c[(c.gis == gis) & (c.status == "ok")].pivot_table(index="cell", columns="survey", values="crest_m_navd88")
    if year not in c.columns or year - 1 not in c.columns:
        return np.nan
    return float((c[year - 1] - c[year]).reindex(d.cell).mean())


# Run: the cells and storms, the sweep, the table, the figure
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--gis", type=int, default=GIS)
    ap.add_argument("--year", type=int, default=YEAR)
    args = ap.parse_args()

    d = G.cells(args.year, (args.gis, args.gis))
    st = E.storms()
    st = st[st.calendar_year == args.year].sort_values("StartTime").reset_index(drop=True)
    rows = [dict(cs=cs, w0_m=w0, mean_lowering_m=lowering(d, st, cs, w0), inside_bounds=True)
            for w0 in W0_VALUES for cs in CS_VALUES]
    rows += [dict(cs=CS_LARSON_LAB, w0_m=w0, mean_lowering_m=lowering(d, st, CS_LARSON_LAB, w0), inside_bounds=False)
             for w0 in W0_VALUES]
    tab = pd.DataFrame(rows)
    tab.to_csv(support_dir(OUT) / "sensitivity.csv", index=False)
    obs = lidar_change(d, args.year, args.gis)
    cal = tab[(tab.cs == CS_CALIBRATED) | np.isclose(tab.cs, CS_CALIBRATED)]
    base = lowering(d, st, CS_CALIBRATED, W0_CALIBRATED)
    lab45 = float(tab[(tab.cs == CS_LARSON_LAB) & (tab.w0_m == W0_CALIBRATED)].mean_lowering_m.iloc[0])
    wide = tab[tab.inside_bounds]
    per_cs = wide[wide.w0_m == W0_CALIBRATED].set_index("cs").mean_lowering_m
    w0_range = [lowering(d, st, CS_CALIBRATED, w) for w in W0_VALUES]
    print(f"GIS {args.gis} {args.year}: calibrated (Cs {CS_CALIBRATED:g}, W0 {W0_CALIBRATED:g}) {base:.2f} m; "
          f"lidar {obs:.2f} m; Cs 0-{DuneCollisionManager.CS_MAX:g} at W0 {W0_CALIBRATED:.0f}: {per_cs.min():.2f}-{per_cs.max():.2f} m; "
          f"W0 {min(W0_VALUES):.0f}-{max(W0_VALUES):.0f} at the calibrated Cs: {min(w0_range):.2f}-{max(w0_range):.2f} m; Larson lab at W0 {W0_CALIBRATED:.0f}: {lab45:.2f} m")

    apply_style()
    fig, (a, b) = plt.subplots(1, 2, figsize=figsize("double", 0.42), layout="constrained")
    for w0, col in zip(W0_VALUES, W0_RAMP):
        s = wide[wide.w0_m == w0]
        a.plot(s.cs * 1e4, s.mean_lowering_m, color=col, lw=1.8, label=f"W0 {w0:.0f} m")
    a.axhline(obs, color=C["REF"], lw=1.0, ls="--")
    a.text(0.05, obs + 0.03, f"lidar {args.year - 1}-{args.year}: {obs:.2f} m", fontsize=7.5, color=C["REF"], va="bottom")
    a.plot(CS_CALIBRATED * 1e4, base, "o", ms=6, color=C["ACCENT"], mec="white", zorder=5)
    a.annotate("hindcast", (CS_CALIBRATED * 1e4, base), xytext=(CS_CALIBRATED * 1e4 + 0.35, base - 0.45),
               fontsize=7.5, arrowprops=dict(arrowstyle="-", lw=0.6))
    a.axvline(CS_LARSON_FIELD * 1e4, color=INK_MUTED, lw=0.7, ls=":")
    a.text(CS_LARSON_FIELD * 1e4 + 0.04, 0.02, "Larson\nfield", fontsize=7, color=INK_MUTED, va="bottom")
    a.set(xlim=(0, DuneCollisionManager.CS_MAX * 1e4), ylim=(0, None), xlabel="Cs (in units of 1e-4)",
          ylabel=f"Mean crest lowering over {args.year} (m)")
    a.legend(loc="upper left", frameon=False, fontsize=7.5)
    open_frame(a)
    _title(a, 0, f"GIS {args.gis}, the {args.year} storms")

    wr = np.linspace(2, 110, 109)
    for cs, col in zip(CS_LINES, ("#c2a5cf", C["ACCENT"], "#40004b", "0.35")):
        b.plot(wr, [lowering(d, st, cs, w) for w in wr], color=col, lw=1.8,
               label=f"Cs {cs * 1e4:.2f}e-4" + (" (hindcast)" if cs == CS_CALIBRATED else
                                       " (tested max)" if cs == DuneCollisionManager.CS_MAX else
                                       " (Larson field)" if cs == CS_LARSON_FIELD else ""))
    b.plot(wr, [lowering(d, st, CS_LARSON_LAB, w) for w in wr], color="0.35", lw=1.2, ls="--",
           label=f"Cs {CS_LARSON_LAB * 1e4:.0f}e-4 (Larson lab, outside bounds)")
    b.axhline(obs, color=C["REF"], lw=1.0, ls="--")
    bw = d.beach_width_m
    b.axvspan(bw.min(), bw.max(), color="0.93", lw=0, zorder=0)
    b.text((bw.min() + bw.max()) / 2, 1.9, f"GIS {args.gis}\nbeaches", ha="center", fontsize=7, color=INK_MUTED, va="bottom")
    b.set(xlim=(0, 110), ylim=(0, None), xlabel="Beach limit W0 (m)", ylabel=f"Mean crest lowering over {args.year} (m)")
    b.legend(loc="upper left", frameon=False, fontsize=7)
    open_frame(b)
    _title(b, 1, "Against the beach limit")
    E.CS, E.W0 = CS_CALIBRATED, W0_CALIBRATED

    caption(fig, f"**Tests:** how much the dune collision module's crest loss depends on its two parameters, Cs and "
            f"the beach limit W0; offline, no model runs. **How to read:** GIS {args.gis}'s {len(d)} dune cells start "
            f"{args.year} at their {d.flown.dropna().iloc[0]} lidar crest, with the hindcast's toes and beach widths "
            f"({bw.min():.0f}-{bw.max():.0f} m from the shoreline to the toe, as dune_toe_cells.csv measures them since "
            f"2026-11-03; shaded in b), and the {len(st)} storms of {args.year} in the hindcast "
            f"storm series lower them through the module's own functions (the year's dune growth left out). (a) the "
            f"mean crest lowering over the year against Cs, from 0 to the module's tested maximum "
            f"{DuneCollisionManager.CS_MAX:g}, for W0 of {', '.join(f'{w:.0f}' for w in W0_VALUES)} m; dot, the "
            f"hindcast's Cs {CS_CALIBRATED:g} and W0 {W0_CALIBRATED:.0f} m; dotted, Larson et al.'s (2004) US East "
            f"Coast value {CS_LARSON_FIELD:g}; dashed green, the lidar's mean crest change on the same cells. (b) the "
            f"same against W0 for four Cs, and dashed grey Larson's laboratory {CS_LARSON_LAB:g}, outside the "
            f"module's bounds. **Shows:** the loss grows almost linearly with Cs until the lowest dunes reach their "
            f"toes: {per_cs.min():.2f}-{per_cs.max():.2f} m over Cs 0-{DuneCollisionManager.CS_MAX:g} at W0 {W0_CALIBRATED:.0f} m; "
            f"the hindcast's setting gives {base:.2f} m against the lidar's {obs:.2f} m. W0 matters as much as Cs on "
            f"these {bw.min():.0f}-{bw.max():.0f} m beaches: {min(w0_range):.2f}-{max(w0_range):.2f} m over W0 {min(W0_VALUES):.0f}-{max(W0_VALUES):.0f} m at the hindcast's Cs. "
            f"Larson's laboratory value would take {lab45:.2f} m, most dunes to their toes.")
    save(fig, OUT / "collision_sensitivity", close=True)


if __name__ == "__main__":
    main()
