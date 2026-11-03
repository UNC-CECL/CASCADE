"""
How the dune collision module (cascade/dune_collision.py) behaves: three figures from its own functions.

    python collision_explained_figures.py
        ->  figures/collision_1_regimes.png
            figures/collision_2_beach.png
            figures/collision_3_in_the_model_year.png

Nothing is fitted or run: one example dune cell (a Buxton dune before Dennis, toe 2.4 m and crest 7.3 m
NAVD88) is put through beach_factor, collision_cells, larson_volume and lower_dunes at the hindcast's
settings (Cs 1.02e-4, beach limit 25 m on the dry beach, beach start to toe), with storms from the 1997-2009 hindcast storm series. Takes a
few seconds. Captions in figures/supporting/CAPTIONS.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-03
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402

_HERE = Path(__file__).resolve()
_HAT_DIR = next(_p for _p in _HERE.parents if (_p / "hatteras_paths.py").exists())
_REPO = next(_p for _p in _HERE.parents if (_p / "pyproject.toml").exists())
for _p in (_HAT_DIR, _REPO):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))
from hatteras_paths import WORKSPACE, workspace_twin  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C, INK, INK_MUTED, _title, apply_style, caption, figsize, open_frame, save)
from cascade.dune_collision import beach_factor, collision_cells, larson_volume, lower_dunes  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
OUT = workspace_twin(_HERE.parent) / "figures"
STORMS = (WORKSPACE / "data" / "hatteras_init" / "3-env-forcings" / "3-storms" / "hindcast_storms"
          / "1997_2009" / "1997_2009_storms_v3_240_summary.csv")
MHW = 0.36                      # m NAVD88
BERM = 1.7                      # m NAVD88, Barrier3D's berm on the Hatteras domains
TOE = 2.4                       # m NAVD88, the example cell's lidar toe (GIS 10-12)
ROWS = (5.6, 4.6)               # m above the berm: front and back dune rows, front crest 7.3 m NAVD88
CS = 1.02e-4                    # Larson et al. (2004) impact coefficient, the hindcast's (1.24e-4 until 2026-11-03)
W0 = 25.0                       # m, dry-beach width (start to toe) at which collision loss stops (45 m to the crest before)
DUNE_W = 20.0                   # m, the two 10 m dune cells
RESTART = 0.075                 # m, Barrier3D's smallest dune height
WIDTH = 10.0                    # m, the example cell's dry beach (GIS 10-12: 0-20 m, median 10)
WIDTHS = (0.0, 8.0, 16.0, 25.0)  # m, dry-beach widths compared in figure 2
YEAR = 1999                     # the storm year of figure 3
SWASH, OVERWASH = 2.0, 7.8      # m NAVD88, example water levels below the toe and above the crest
C_MOD, C_MOD_FILL = C["ACCENT"], C["ACCENT_FILL"]     # the module: what is under test
C_B3D = C["BASE"]                                     # Barrier3D's own behaviour, unchanged
C_SAND, C_WATER = "#e8dcc0", C["WATER"]
# -----------------------------------------------------------------------------


# m NAVD88 to the module's dam MHW, and back
def dam(z):
    return (np.asarray(z, dtype=float) - MHW) / 10.0


def m(z_dam):
    return np.asarray(z_dam, dtype=float) * 10.0 + MHW


# One storm on one or many identical cells: the crest after it (m NAVD88), as DuneCollisionManager.update
def crest_after(r_high, period, duration, width, rows=ROWS, crest_prestorm=None):
    r = np.atleast_1d(np.asarray(r_high, dtype=float))
    n = len(r)
    rows_dam = np.tile(np.asarray(rows) / 10.0, (n, 1))
    crest = rows_dam.max(axis=1) + dam(BERM)
    pre = crest if crest_prestorm is None else np.full(n, dam(crest_prestorm))
    toe = np.full(n, dam(TOE))
    beach = beach_factor(np.full(n, width / 10.0), W0 / 10.0)
    out = m(crest)
    for i in range(n):
        if collision_cells(pre[i:i + 1], crest[i:i + 1], toe[i:i + 1], dam(r[i]), beach[i:i + 1]).any():
            v = larson_volume(dam(r[i]), toe[i:i + 1], period, duration, CS, beach[i:i + 1])
            new = lower_dunes(rows_dam[i:i + 1], v, toe[i:i + 1], dam(BERM), DUNE_W / 10.0, RESTART / 10.0)
            out[i] = m(new.max() + dam(BERM))
    return out


# The rows after one storm (m above the berm), for the profile sketch
def rows_after(r_high, period, duration, width):
    rows_dam = np.asarray(ROWS)[None, :] / 10.0
    toe = np.array([dam(TOE)])
    beach = beach_factor(np.array([width / 10.0]), W0 / 10.0)
    v = larson_volume(dam(r_high), toe, period, duration, CS, beach)
    return lower_dunes(rows_dam, v, toe, dam(BERM), DUNE_W / 10.0, RESTART / 10.0)[0] * 10.0


# The storm series' storms, with R_high in m NAVD88
def storms():
    s = pd.read_csv(STORMS, parse_dates=["StartTime"])
    return s.assign(r_m=s.Rhigh * 10.0 + MHW)


# A schematic cross-shore profile: interior, two dune rows, toe, beach to MHW (x seaward, m)
def profile(rows, width):
    front, back = BERM + rows[0], BERM + rows[1]
    toe_x = 40.0
    x = [-10, 0, 4, 10, 16, 22, 28, 34, toe_x, toe_x + width, toe_x + width + 25]
    z = [BERM, BERM, back * 0.8 + BERM * 0.2, back, (back + front) / 2, front, front, (front + TOE) / 2,
         TOE, MHW, MHW - 1.5]
    return np.array(x), np.array(z), toe_x


# Figure 1: the three regimes on a profile, and the crest after one storm against its water level
def fig_regimes(st):
    dennis = st[(st.calendar_year == 1999) & (st.r_m == st[st.calendar_year == 1999].r_m.max())].iloc[0]
    T, D = dennis.period, dennis.duration
    crest0 = BERM + max(ROWS)

    fig, (a, b) = plt.subplots(1, 2, figsize=figsize("double", 0.42), layout="constrained",
                               gridspec_kw=dict(width_ratios=(1.25, 1)))
    # (a) the profile, the three water levels, the lowered dune
    x, z, toe_x = profile(ROWS, 0.0 + 25.0)
    a.fill_between(x, z, MHW - 2, color=C_SAND, lw=0)
    after = rows_after(dennis.r_m, T, D, 0.0)
    xa, za, _ = profile(after, 25.0)
    a.fill_between(x, z, np.minimum(z, za), where=z > za + 1e-6, color=C_MOD_FILL, lw=0,
                   label=f"lost to collision in one storm ({dennis.r_m:.2f} m, {D:.0f} h)")
    a.plot(x, z, color=INK, lw=0.9)
    a.plot(xa, za, color=C_MOD, lw=1.2, ls="--")
    for zl, name, col in ((SWASH, "swash: below the toe, nothing", INK_MUTED),
                          (dennis.r_m, "collision: toe to crest, this module", C_MOD),
                          (OVERWASH, "overwash: above the crest, Barrier3D", C_B3D)):
        # the water meets the seaward face or the beach; above the crest it runs across
        xl = np.interp(zl, z[6:][::-1], x[6:][::-1]) if zl < crest0 else -10
        a.plot([xl, 92], [zl, zl], color=col, lw=1.4)
        a.text(92, zl + 0.12, name, ha="right", va="bottom", fontsize=7.5, color=col)
    a.plot(toe_x, TOE, "o", ms=4, color=INK)
    a.annotate("toe", (toe_x, TOE), xytext=(toe_x + 4, TOE - 0.9), fontsize=7.5, color=INK,
               arrowprops=dict(arrowstyle="-", lw=0.6, color=INK))
    a.annotate("crest", (22, crest0), xytext=(4, crest0 + 0.7), fontsize=7.5, color=INK,
               arrowprops=dict(arrowstyle="-", lw=0.6, color=INK))
    a.text(toe_x + 12, MHW + 0.25, "beach", fontsize=7.5, color=INK_MUTED)
    a.set(xlim=(-10, 92), ylim=(MHW - 1.5, 9.4), xlabel="Cross-shore distance (m, seaward →)",
          ylabel="Elevation (m NAVD88)")
    a.legend(loc="upper left", frameon=False, fontsize=7.5)
    open_frame(a)
    _title(a, 0, "The three regimes")

    # (b) crest after one Dennis-length storm against its water level, no beach
    r = np.linspace(1.5, 8.5, 400)
    ca = crest_after(r, T, D, 0.0)
    ca[r > crest0] = np.nan
    b.axvspan(1.5, TOE, color="0.95", lw=0)
    b.axvspan(crest0, 8.5, color="0.90", lw=0)
    b.plot(r, ca, color=C_MOD, lw=2)
    b.axhline(TOE, color=INK_MUTED, lw=0.8, ls=":")
    b.text(8.4, TOE + 0.1, "floor: the toe", ha="right", va="bottom", fontsize=7.5, color=INK_MUTED)
    b.plot(dennis.r_m, crest_after([dennis.r_m], T, D, 0.0)[0], "o", ms=5, color=C_MOD, mec="white")
    b.annotate(f"Dennis {dennis.r_m:.2f} m", (dennis.r_m, crest_after([dennis.r_m], T, D, 0.0)[0]),
               xytext=(dennis.r_m + 0.35, 7.6), fontsize=7.5, arrowprops=dict(arrowstyle="-", lw=0.6))
    b.text((1.5 + TOE) / 2, 1.0, "swash", ha="center", fontsize=7.5, color=INK_MUTED)
    b.text((crest0 + 8.5) / 2, 1.0, "overwash\n(Barrier3D)", ha="center", fontsize=7.5, color=C_B3D)
    b.text((TOE + crest0) / 2, 1.0, "collision", ha="center", fontsize=7.5, color=C_MOD)
    b.set(xlim=(1.5, 8.5), ylim=(0.6, 8.2), xlabel="Storm's peak water level R_high (m NAVD88)",
          ylabel="Crest after the storm (m NAVD88)")
    open_frame(b)
    _title(b, 1, "One storm, one dune, no beach")

    low = crest0 - crest_after([dennis.r_m], T, D, 0.0)[0]
    floor_r = r[np.nanargmin(np.abs(ca - TOE - 1e-3))] if np.nanmin(ca) <= TOE + 1e-3 else np.nan
    caption(fig, f"**Tests:** nothing new; this figure explains which storms the dune collision module "
            f"(`cascade/dune_collision.py`, `DuneCollisionManager`) acts on, and how far it lowers a dune. "
            f"**How to read:** (a) a schematic cross-shore profile of one Barrier3D dune cell, ocean to the "
            f"right: two 10 m dune rows {ROWS[0]:.1f} and {ROWS[1]:.1f} m above the {BERM} m berm (crest "
            f"{crest0:.1f} m NAVD88, as GIS 12 before Dennis in the hindcast), toe {TOE} m NAVD88 (the lidar "
            f"toe at GIS 10-12). Three water levels (Sallenger 2000): below the toe nothing happens; between "
            f"toe and crest is collision, the module's regime; above the pre-storm crest Barrier3D's own "
            f"overwash lowers the dune and the module leaves the cell alone. Purple: the sand one storm the "
            f"size of Dennis in the storm series ({dennis.r_m:.2f} m held for {D:.0f} h, T {T:.1f} s) removes "
            f"with no beach in front; both rows are lowered in proportion. (b) the crest after that storm "
            f"against its peak water level R_high, the same cell, Cs {CS:g}, no beach. **Shows:** the loss "
            f"starts at the toe and grows with the square of the water above it (Larson et al. 2004, "
            f"dV = 4 Cs (R_high - toe)^2 t / T, spread over the 20 m dune): Dennis lowers this crest "
            f"{low:.2f} m; from about {floor_r:.1f} m the storm takes the dune to its toe, the floor. Above "
            f"the crest the curve stops: that cell is overwashed, and Barrier3D lowers it as before.")
    save(fig, OUT / "collision_1_regimes", close=True)
    return dennis, low


# Figure 2: the beach in front of the dune scales the loss down to nothing at W0
def fig_beach(dennis):
    T, D = dennis.period, dennis.duration
    crest0 = BERM + max(ROWS)
    fig, (a, b) = plt.subplots(1, 2, figsize=figsize("double", 0.40), layout="constrained")

    w = np.linspace(0, 60, 200)
    a.plot(w, beach_factor(w / 10.0, W0 / 10.0), color=C_MOD, lw=2)
    a.axvline(W0, color=INK_MUTED, lw=0.8, ls=":")
    a.text(W0 + 1, 0.85, f"W0 = {W0:.0f} m:\nno collision loss\nbehind a wider beach", fontsize=7.5,
           color=INK_MUTED, va="top")
    ramp = ("#c2a5cf", "#9970ab", "#7b3294", "#40004b")
    for wi, col in zip(WIDTHS, ramp):
        a.plot(wi, beach_factor(np.array([wi / 10.0]), W0 / 10.0)[0], "o", ms=6, color=col, mec="white")
    a.set(xlim=(0, 60), ylim=(-0.03, 1.05), xlabel="Dry beach width W, beach start to toe (m)",
          ylabel="Share of Larson's volume lost\nmax(0, 1 − W / W0)")
    open_frame(a)
    _title(a, 0, "The beach factor")

    r = np.linspace(TOE - 0.3, crest0, 300, endpoint=False)   # at the crest it is overwash
    lows = {}
    for wi, col in zip(WIDTHS, ramp):
        ca = crest_after(r, T, D, wi)
        b.plot(r, ca, color=col, lw=2, label=f"{wi:.0f} m beach")
        lows[wi] = crest0 - crest_after([dennis.r_m], T, D, wi)[0]
    b.axhline(TOE, color=INK_MUTED, lw=0.8, ls=":")
    b.axvline(dennis.r_m, color=INK_MUTED, lw=0.8)
    b.text(dennis.r_m + 0.08, 7.75, f"Dennis {dennis.r_m:.2f} m", fontsize=7.5, color=INK_MUTED)
    b.set(xlim=(TOE - 0.3, crest0), ylim=(TOE - 0.4, 8.0), xlabel="Storm's peak water level R_high (m NAVD88)",
          ylabel="Crest after the storm (m NAVD88)")
    b.legend(loc="center left", bbox_to_anchor=(0.0, 0.42), frameon=False)
    open_frame(b)
    _title(b, 1, "The same storm behind four beaches")

    caption(fig, f"**Tests:** nothing new; this figure explains the beach factor of the dune collision module, "
            f"the part that makes it fit the Hatteras lidar (`../1-observations/`: crest loss in collision "
            f"follows beach width, not R_high/D_high). **How to read:** (a) the share of Larson's volume a "
            f"dune loses against the width of the beach in front of it, max(0, 1 - W / W0) with the "
            f"hindcast's W0 = {W0:.0f} m; dots, the four widths of (b). (b) the crest of the figure-1 cell "
            f"(crest {crest0:.1f} m, toe {TOE} m NAVD88) after one Dennis-length storm ({D:.0f} h, T {T:.1f} s) "
            f"against its peak water level, behind beaches of {', '.join(f'{x:.0f}' for x in WIDTHS)} m; "
            f"each W is held for the run (the cell's median over the surveys). **Shows:** the beach does not "
            f"change which storms reach the dune, only how much they take: Dennis lowers this crest "
            + ", ".join(f"{lows[x]:.2f} m behind {x:.0f} m" for x in WIDTHS)
            + f". Behind a dry beach of {W0:.0f} m or more the module does nothing, so the tall dunes behind wide beaches "
            "keep their height, as the lidar shows.")
    save(fig, OUT / "collision_2_beach", close=True)


# Figure 3: where the module acts in Barrier3D's year, and one year of storms on the example cell
def fig_year(st):
    crest0 = BERM + max(ROWS)
    fig, (a, b) = plt.subplots(1, 2, figsize=figsize("double", 0.46), layout="constrained",
                               gridspec_kw=dict(width_ratios=(1.05, 1.15)))
    # (a) the order within one Barrier3D year
    a.set(xlim=(-0.12, 1.08), ylim=(0, 1.06)); a.axis("off")
    steps = [(0.90, "Dune growth (Barrier3D)", C_B3D, "then the module keeps each\ncell's pre-storm crest"),
             (0.68, "Storm k: overwash (Barrier3D)", C_B3D, "lowers cells with\nR_high > pre-storm crest"),
             (0.46, "Storm k: collision (module)", C_MOD, "lowers cells with\ntoe < R_high < crest"),
             (0.24, "Dune diffusion (Barrier3D)", C_B3D, "then storm k + 1,\non the lowered dunes"),
             (0.04, "Rest of the year", C_B3D, "interior and shoreline,\nthen BRIE alongshore")]
    for y, text, col, note in steps:
        a.add_patch(FancyBboxPatch((0.0, y), 0.62, 0.13, boxstyle="round,pad=0.01,rounding_size=0.02",
                                   fc="white", ec=col, lw=1.6 if col == C_MOD else 0.9))
        a.text(0.31, y + 0.065, text, ha="center", va="center", fontsize=7,
               color=col if col == C_MOD else INK, fontweight="bold" if col == C_MOD else "normal")
        if note:
            a.text(0.65, y + 0.065, note, ha="left", va="center", fontsize=7, color=INK_MUTED)
    for y0, y1 in ((0.90, 0.81), (0.68, 0.59), (0.46, 0.37), (0.24, 0.17)):
        a.annotate("", (0.31, y1), (0.31, y0), arrowprops=dict(arrowstyle="-|>", lw=0.8, color=INK))
    # the storm loop: back from dune diffusion to the next storm's overwash
    a.plot([-0.005, -0.05, -0.05], [0.305, 0.305, 0.745], color=INK_MUTED, lw=0.9)
    a.annotate("", (-0.005, 0.745), (-0.05, 0.745), arrowprops=dict(arrowstyle="-|>", lw=0.9, color=INK_MUTED))
    a.text(-0.07, 0.525, "next storm", ha="right", va="center", fontsize=7, color=INK_MUTED, rotation=90)
    _title(a, 0, "Where the module acts")

    # (b) the year's storms on the example cell
    yr = st[st.calendar_year == YEAR].sort_values("StartTime").reset_index(drop=True)
    crests, c = [crest0], crest0
    rows = np.array(ROWS, dtype=float)
    for _, s in yr.iterrows():
        new = crest_after([s.r_m], s.period, s.duration, WIDTH, rows=tuple(rows), crest_prestorm=crest0)[0]
        if new < c:
            rows = rows * (new - BERM) / (c - BERM)
        c = new
        crests.append(c)
    t = yr.StartTime
    t0 = pd.Timestamp(f"{YEAR}-01-01")
    t_end = pd.Timestamp(f"{YEAR}-12-31")
    reach = yr.r_m > TOE
    b.vlines(t[~reach], MHW, yr.r_m[~reach], color=C_WATER, lw=1.6)
    b.vlines(t[reach], MHW, yr.r_m[reach], color=C["LATE"], lw=1.6)
    b.step(pd.concat([pd.Series([t0]), t, pd.Series([t_end])]), crests + [crests[-1]], where="post",
           color=C_MOD, lw=2, label=f"crest, {WIDTH:.0f} m beach")
    b.axhline(TOE, color=INK_MUTED, lw=0.8, ls=":")
    b.text(t_end, TOE + 0.08, "toe", ha="right", va="bottom", fontsize=7.5, color=INK_MUTED)
    big = yr.r_m.idxmax()
    b.annotate(f"Dennis {yr.r_m[big]:.2f} m, {yr.duration[big]:.0f} h: −{crests[big] - crests[big + 1]:.2f} m",
               (t[big], crests[big + 1]), xytext=(pd.Timestamp(f"{YEAR}-01-20"), 5.4), fontsize=7.5,
               arrowprops=dict(arrowstyle="-", lw=0.6))
    b.plot([], [], color=C["LATE"], lw=1.6, label="storm R_high, reaches the toe")
    b.plot([], [], color=C_WATER, lw=1.6, label="storm R_high, below the toe")
    b.set(xlim=(t0, t_end), ylim=(0, 8.2), ylabel="Elevation (m NAVD88)")
    b.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%b"))
    b.legend(loc="center right", bbox_to_anchor=(1.0, 0.56), frameon=False, fontsize=7.5)
    open_frame(b)
    _title(b, 1, f"The {YEAR} storms on one dune")

    acted = int(sum(np.diff(crests) < 0))
    caption(fig, f"**Tests:** nothing new; this figure shows where the dune collision module acts inside a "
            f"Barrier3D year, and what a year of storms does to one dune. **How to read:** (a) the order "
            f"within one Barrier3D year: the module stands in for two Barrier3D methods on its domain without "
            f"changing them. After Barrier3D's dune growth it keeps each cell's pre-storm crest; after each "
            f"storm's overwash lowering, and before Barrier3D's dune diffusion, it lowers the cells that storm "
            f"left in collision. Purple, the module; grey, Barrier3D unchanged. (b) the {YEAR} storms of the "
            f"hindcast storm series ({len(yr)} storms; stems, each storm's peak R_high) on the figure-1 cell "
            f"(crest {crest0:.1f} m, toe {TOE} m NAVD88) behind a {WIDTH:.0f} m dry beach (GIS 10-12 have 0-20 m), "
            f"Cs {CS:g}, W0 {W0:.0f} m; the purple step is the crest after each storm. **Shows:** "
            f"{int(reach.sum())} of the {len(yr)} storms reach the toe and {acted} lower the dune; Dennis does "
            f"most of it, and the storms after it see the lowered dune. Over the year the crest falls from "
            f"{crest0:.2f} to {crests[-1]:.2f} m. The module only lowers: regrowth is Barrier3D's dune growth, "
            f"once a year before the storms.")
    save(fig, OUT / "collision_3_in_the_model_year", close=True)
    return crests


# Run: style, the three figures
def main() -> None:
    apply_style()
    st = storms()
    dennis, low = fig_regimes(st)
    fig_beach(dennis)
    crests = fig_year(st)
    print(f"Dennis on the example cell, no beach: -{low:.2f} m; {YEAR} crest {crests[0]:.2f} -> {crests[-1]:.2f} m")
    print(f"figures in {OUT}")


if __name__ == "__main__":
    main()
