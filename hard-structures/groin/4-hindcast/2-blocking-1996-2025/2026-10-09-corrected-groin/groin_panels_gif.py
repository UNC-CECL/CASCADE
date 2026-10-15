#!/usr/bin/env python3
"""
The groin's and the fills' effect on the modelled shoreline, year by year, as stacked plan views.

    python groin_panels_gif.py   ->  figures/corrected_groin_panels_<calibration|test>.gif

One panel per run, stacked, each a plan view along GIS 2-12: the shoreline's cross-shore position in the
model's frame, rotated so the coast's general trend (a straight-line fit to the start line) runs left to
right, land up and ocean down, so accretion moves the line down and erosion up. Stacking gives each panel
the full width, which keeps the cross-shore exaggeration to about 4-5 times (printed on each GIF, with
scale bars). Drawn: the start line (observed and model), the observed end line (the target), and this
year's modelled line, with accretion (yellow) and erosion (red) shaded between the start and the model.
Test, 2009-2025: no groin and no fills; fills only; fills and the corrected groin (0.10, off from 2017).
Calibration, 1996-2009, has no fills: no groin; the corrected groin (0.16). Calendar-year setup: row k of
a run is calendar year start + k. GIS 1 is left out (the solved end moves it 150 m). `render` is reused by
../2026-10-13-pregroin-planform/.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-13
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.animation import PillowWriter  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

import groin_on_off_gif as G  # noqa: E402
from site_layer.hat_figure_style import C, INK, INK_MUTED, apply_style, figsize, open_frame, panel_title  # noqa: E402
from site_layer.hatteras_site_config import last_model_year  # noqa: E402

P, rf, pv, SR, cgm = G.P, G.rf, G.pv, G.SR, G.cgm

# --- CONFIG ------------------------------------------------------------------
GIS = np.arange(2, 13)                     # domains drawn, south to north
DOMAIN_M = 500.0                           # alongshore length of one domain
ACC_C, ERO_C = "#f2c12e", C["EARLY"]       # accretion yellow, erosion red
LAND_C, SEA_C = "#f4f1ea", "#e3eef6"       # the start's land and ocean, kept pale so the change stands out
OBS_C = C["LOCATOR"]                       # observations: the house teal
GROIN_X, GROIN_LEN = 5.5, 35.0             # the groin face between GIS 5 and 6; drawn length (m, seaward)
FILLS = {2017: "2017", 2022: "2022"}       # Buxton fills, GIS 6-16
PANEL_H = 1.25                             # inches per stacked panel
FPS = 2
HOLD = 6
# -----------------------------------------------------------------------------


# Take the start line's straight-line trend off every line, so the coast's general direction runs left to right
def rotate(lines, start):
    trend = np.polyval(np.polyfit(GIS, start, 1), GIS)
    return [np.asarray(v) - trend for v in lines]


# The stacked GIF. runs: (title, [year x domain] positions after rotate, has fills, groin on in frame k or None)
def render(runs, x0, target, start_label, target_label, frames, title, status, out, fills_from=None):
    allv = np.concatenate([m.ravel() for _, m, _, _ in runs] + [target, x0])
    lo, hi = np.nanmin(allv) - 35, np.nanmax(allv) + 25
    n = len(runs)
    fig, axes = plt.subplots(n, 1, figsize=(figsize("double")[0], 1.5 + PANEL_H * n + 0.75), sharex=True)
    fig.subplots_adjust(left=0.1, right=0.98, top=1 - 0.8 / (1.5 + PANEL_H * n + 0.75),
                        bottom=0.95 / (1.5 + PANEL_H * n + 0.75), hspace=0.32)
    axes = np.atleast_1d(axes)
    h = [Line2D([], [], color=INK_MUTED, lw=1.0, ls=(0, (3, 2)), label=start_label),
         Line2D([], [], color=OBS_C, lw=1.2, ls="--", marker="s", ms=3.5, mfc="white", label=target_label),
         Line2D([], [], color=INK, lw=1.8, label="modelled shoreline"),
         Patch(facecolor=ACC_C, label="accretion since the start"),
         Patch(facecolor=ERO_C, alpha=0.8, label="erosion since the start"),
         Line2D([], [], color=INK, lw=3.0, solid_capstyle="butt", label="Buxton groins")]
    fig.legend(handles=h, loc="lower center", ncol=3, fontsize=6.6, frameon=False)
    fig.text(0.1, 0.995, title, fontsize=8.5, ha="left", va="top", color=INK)
    year_t = fig.text(0.98, 0.995, "", fontsize=15, fontweight="bold", ha="right", va="top", color=INK)
    state_t = fig.text(0.1, 0.995 - 0.3 / fig.get_figheight(), "", fontsize=7, ha="left", va="top", color=INK_MUTED)
    # Exaggeration of the cross-shore axis against the alongshore one, from the drawn axes
    w, hgt = axes[0].get_position().width * fig.get_figwidth(), axes[0].get_position().height * fig.get_figheight()
    exag = ((len(GIS)) * DOMAIN_M / w) / ((hi - lo) / hgt)

    def draw(k):
        for i, (ax, (ttl, m, fills, groin_on)) in enumerate(zip(axes, runs)):
            ax.clear()
            ax.fill_between(GIS, lo, x0, color=SEA_C, lw=0, zorder=0)
            ax.fill_between(GIS, x0, hi, color=LAND_C, lw=0, zorder=0)
            ax.fill_between(GIS, x0, m[k], where=m[k] < x0, color=ACC_C, lw=0, interpolate=True, zorder=1)
            ax.fill_between(GIS, x0, m[k], where=m[k] > x0, color=ERO_C, alpha=0.8, lw=0, interpolate=True, zorder=1)
            ax.plot(GIS, x0, color=INK_MUTED, lw=1.0, ls=(0, (3, 2)), zorder=3)
            ax.plot(GIS, target, "--", color=OBS_C, lw=1.2, marker="s", ms=3.5, mfc="white", zorder=4)
            ax.plot(GIS, m[k], "-", color=INK, lw=1.8, zorder=5)
            g0 = np.interp(GROIN_X, GIS, x0)
            on = groin_on(k) if groin_on else None
            if on is None:
                ax.axvline(GROIN_X, color=INK_MUTED, lw=0.6, ls=(0, (1, 2)), zorder=2)
                ax.text(GROIN_X + 0.1, hi - 4, "groin site", fontsize=6.2, color=INK_MUTED, va="top")
            else:                                             # the structure, from the start shoreline seaward
                ax.plot([GROIN_X, GROIN_X], [g0, g0 - GROIN_LEN], color=INK if on else INK_MUTED, lw=3.0,
                        solid_capstyle="butt", alpha=1 if on else 0.5, zorder=6)
                ax.text(GROIN_X - 0.12, g0 - GROIN_LEN, "groin" if on else "groin (off)", fontsize=6.2,
                        color=INK if on else INK_MUTED, ha="right", va="center", zorder=6)
            placed = [nm for fy, nm in FILLS.items() if fills and fills_from is not None and fills_from(k) >= fy]
            if placed:
                ax.plot([GROIN_X + 0.5, GIS[-1] + 0.45], [hi - 6, hi - 6], color=INK, lw=1.0, zorder=6)
                ax.text((GROIN_X + 0.5 + GIS[-1] + 0.45) / 2, hi - 9, f"{' and '.join(placed)} fill"
                        f"{'s' if len(placed) > 1 else ''} (GIS 6-16)", fontsize=6.2, ha="center", va="top",
                        color=INK, zorder=6)
            if i < n - 1:
                ax.text(GIS[0] - 0.4, lo + 5, "ocean", fontsize=6.5, color=INK_MUTED, style="italic", va="bottom")
            ax.text(GIS[0] - 0.4, hi - 4, "land", fontsize=6.5, color=INK_MUTED, style="italic", va="top")
            ax.set_xlim(GIS[0] - 0.5, GIS[-1] + 0.5)
            ax.set_ylim(lo, hi)
            ax.set_xticks(GIS)
            ax.tick_params(labelsize=6.8)
            open_frame(ax)
            ax.set_title(panel_title(chr(ord("a") + i), ttl), fontsize=7.8, loc="left", color=INK, pad=3)
        axes[-1].set_xlabel("GIS domain (500 m each), south to north; drift runs north to south", fontsize=7.5)
        fig.supylabel("Cross-shore position about the coast's general trend (m; land up)", fontsize=7, x=0.015)
        ax = axes[-1]                                         # scale bars, bottom right
        x1, yb = GIS[-1] + 0.05, lo + 8
        ax.plot([x1 - 1, x1], [yb, yb], color=INK, lw=1.2, zorder=7)
        ax.text(x1 - 0.5, yb + 2, "500 m", fontsize=6, ha="center", va="bottom", color=INK, zorder=7)
        ax.plot([x1, x1], [yb, yb + 50], color=INK, lw=1.2, zorder=7)
        ax.text(x1 + 0.06, yb + 25, "50 m", fontsize=6, ha="left", va="center", color=INK, zorder=7)
        ax.text(GIS[0] - 0.4, lo + 5, f"ocean    (cross-shore exaggerated x{exag:.0f})", fontsize=6.5,
                color=INK_MUTED, style="italic", va="bottom")
        year_t.set_text(frames[k])
        state_t.set_text(status(k))

    writer = PillowWriter(fps=FPS)
    with writer.saving(fig, str(out), dpi=150):
        for k in range(len(frames)):
            draw(k)
            writer.grab_frame()
        for _ in range(HOLD):
            writer.grab_frame()
    plt.close(fig)
    print("wrote", out.relative_to(cgm.REPO), f"(cross-shore x{exag:.1f})")


# The runs drawn in each period: (panel title, shoreline file, has fills, groin strength or None)
def panels(period):
    base = SR.baseline_file(period)
    if period == rf.CAL:
        return [("No groin", base, False, None),
                ("Corrected groin (strength 0.16)", rf.two_number_file(period, 0.16), False, 0.16)]
    nofill = base.parent.parent / base.parent.name.replace("_nourish_", "_nonourish_")
    nofill = nofill / f"{nofill.name}_shoreline_matrix.npy"
    return [("No groin, no fills", nofill, False, None),
            ("Fills, no groin", base, True, None),
            ("Fills and groin (0.10, off from 2017)", rf.two_number_file(period, 0.10), True, 0.10)]


def make(period, data):
    name = "calibration" if period == rf.CAL else "test"
    cgm.observed_change(period)                               # sets the calendar start window
    last = last_model_year(period)
    years = list(range(period, last + 1))
    start = P.start_positions(data, period, data.dom[data.dom.isin(GIS)].index).dropna()
    obs_end = G.observed_profile(data, last, start).reindex(GIS).to_numpy()
    spec = panels(period)
    idx = GIS + P.PAD - 1
    # Positions in the model's cross-shore frame (landward +), so the coast keeps its shape: land up, ocean down
    model = [np.load(path)[: len(years)][:, idx] for _, path, _, _ in spec]
    x0_raw = model[0][0]                                      # every run starts on the same observed line
    *model, x0, target = rotate(model + [x0_raw, x0_raw - obs_end], x0_raw)   # observed change is seaward +
    runs = [(t, m, fills, (None if s is None else (lambda k: period == rf.CAL or years[k] < 2017)))
            for (t, _, fills, s), m in zip(spec, model)]

    def status(k):
        y = years[k]
        if period == rf.CAL:
            return ""
        return "; ".join(["groin on" if y < 2017 else "groin off from the 2017 fill"]
                         + [f"{nm} fill placed" for fy, nm in FILLS.items() if y >= fy])

    render(runs, x0, target, f"{period} shoreline (start)", f"{last} observed shoreline (target)",
           [str(y) for y in years],
           f"Modelled shoreline with and without the groin{' and fills' if period == rf.TEST else ''}, "
           f"{name} run {period}-{last}", status, rf.FIGS / f"corrected_groin_panels_{name}.gif",
           fills_from=lambda k: years[k])


def main():
    apply_style()
    data = pv.Data()
    for period in (rf.CAL, rf.TEST):
        make(period, data)


if __name__ == "__main__":
    main()
