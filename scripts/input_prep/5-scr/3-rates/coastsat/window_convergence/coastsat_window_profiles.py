"""
What does each nested window's alongshore rate profile look like against 1996-2024?

    python scripts/input_prep/5-scr/3-rates/coastsat/window_convergence/coastsat_window_profiles.py
    python scripts/input_prep/5-scr/3-rates/coastsat/window_convergence/coastsat_window_profiles.py --direction forward

Fits every nested window at every transect and draws the profiles; writes the
transects table that coastsat_window_r_bias_rmse.py reads, figures and a README
per direction. Details: scripts/input_prep/5-scr/3-rates/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-02
"""

import argparse
import datetime
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Rule 5: find the root by searching upward.
_REPO = next(_p for _p in Path(__file__).resolve().parents
             if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(_REPO / "scripts"))
from site_layer import hat_observed_rates as obs      # noqa: E402
from site_layer import hat_figure_style as fs         # noqa: E402

# The sibling sweep, for its loader and its per-transect fit
sys.path.insert(0, str(Path(__file__).resolve().parent))
import coastsat_window_convergence as wc              # noqa: E402

import matplotlib                                     # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                       # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, Normalize  # noqa: E402
from matplotlib.cm import ScalarMappable              # noqa: E402


# --- CONFIG ------------------------------------------------------------------
REF_START, REF_END = wc.REF_START, wc.REF_END

# TWO YEARS, not the sweep's five (Hannah, 2026-09-29)
MIN_WINDOW_YEARS = 2

# Windows light (short) to dark (long) on the shoreline-blue ramp
WINDOW_CMAP = LinearSegmentedColormap.from_list(
    "window_length", ["#c6dbef", "#6baed6", fs.C_1997, "#08306b"])
REF_COLOUR = fs.C["ACCENT"]

# THE PANEL FIGURE IS ZOOMED (Hannah, 2026-09-29
PANEL_Y_HALF = 5.0

# The overlay is zoomed too (Hannah, 2026-10-01): the full axis hid the patterns
OVERLAY_Y_HALF = 10.0

# No window is singled out in the panels (Hannah, 2026-09-29)
PANEL_TITLE_PT = 9.0
PANEL_TICK_PT = 8.0
PANEL_LABEL_PT = 9.0

# v2 of the forward panels outlines the chosen window in yellow (Hannah, 2026-10-02)
CHOSEN_WINDOW = {"forward": (1996, 2015), "backward": (2010, 2026)}
CHOSEN_COLOUR = "#ffd400"
CHOSEN_LW = 3.0
# v3: only the windows of these lengths in years of record, the chosen one still outlined (Hannah, 2026-10-02)
V3_YEARS = {"forward": (14, 25), "backward": (14, 25)}
# -----------------------------------------------------------------------------


# The nested family, SHORTEST FIRST, so the reference is always last
def windows_for(direction):
    if direction == "forward":
        return [(REF_START, y, y)
                for y in range(REF_START + MIN_WINDOW_YEARS - 1, REF_END + 1)]
    if direction == "backward":
        return [(y, REF_END, y)
                for y in range(REF_END - MIN_WINDOW_YEARS + 1, REF_START - 1, -1)]
    raise ValueError("direction is 'forward' or 'backward', not {0!r}".format(direction))


# A window as 'start-end'
def window_label(start, end):
    return "{0}–{1}".format(start, end)


# Figure title: the direction and what is drawn (Hannah, 2026-10-01)
def figure_title(direction):
    if direction == "forward":
        head = "Forward from {0}: start fixed at {0}, end moves later".format(REF_START)
    else:
        head = "Backward from {0}: end fixed at {0}, start moves earlier".format(REF_END)
    return "{0}\nShoreline change rate for each window (blue) vs. {1} (purple)".format(
        head, window_label(REF_START, REF_END))


# The sweep

# Every window's domain profile for a direction
def run(direction):
    windows = windows_for(direction)
    picks = wc.all_transects()
    print("  {0}: {1} transects x {2} windows ({3} ... {4})".format(
        direction, len(picks), len(windows),
        window_label(*windows[0][:2]), window_label(*windows[-1][:2])))
    rows = []
    for i, (domain, tid) in enumerate(picks, start=1):
        rows.extend(wc.sweep_one(tid, domain, windows))
        if i % 150 == 0 or i == len(picks):
            print("    {0}/{1} transects".format(i, len(picks)))
    sweep = pd.DataFrame(rows)
    sweep = sweep.merge(wc.alongshore_x(picks)[["transect_id", "x_domain"]],
                        on="transect_id", how="left")
    ref = (sweep[(sweep["start_year"] == REF_START) & (sweep["end_year"] == REF_END)]
           .set_index("transect_id")["lrr_m_yr"])
    sweep["ref_lrr_m_yr"] = sweep["transect_id"].map(ref)
    sweep["diff_m_yr"] = sweep["lrr_m_yr"] - sweep["ref_lrr_m_yr"]
    return windows, sweep


# One row per window
def correlations(sweep, windows):
    out = []
    for start, end, moving in windows:
        w = sweep[(sweep["start_year"] == start) & (sweep["end_year"] == end)]
        ok = w[["lrr_m_yr", "ref_lrr_m_yr"]].dropna()
        r = (float(np.corrcoef(ok["lrr_m_yr"], ok["ref_lrr_m_yr"])[0, 1])
             if len(ok) > 2 else np.nan)
        out.append({
            "window": "{0}_{1}".format(start, end),
            "start_year": start, "end_year": end, "moving_year": moving,
            "n_years": end - start + 1,
            "r_vs_reference": round(r, 4) if r == r else np.nan,
            "n_transects": len(ok),
            "n_transects_no_fit": int(w["lrr_m_yr"].isna().sum()),
            "median_lrr_m_yr": round(float(w["lrr_m_yr"].median()), 4),
            "median_diff_m_yr": round(float(w["diff_m_yr"].median()), 4),
        })
    return pd.DataFrame(out)


# Figures

# x and rate for one window in alongshore order, NaN where there is no fit so the line breaks
def _profile(frame):
    frame = frame.sort_values("x_domain")
    return frame["x_domain"].to_numpy(), frame["lrr_m_yr"].to_numpy()


# Zero line, limits and grid for a profile panel
def _finish_axis(ax, label_towns):
    ax.axhline(0.0, color=fs.C["INK_MUTED"], lw=0.5, zorder=1)
    ax.set_xlim(0.5, 90.5)
    ax.grid(True, axis="y", alpha=0.6)
    fs.town_bands(ax, label=label_towns)


# Every window's profile over the reference
def draw_overlay(sweep, windows, direction, out_dir):
    fs.apply_style()
    fig, ax = plt.subplots(figsize=fs.figsize("double", height=4.4), layout="constrained")
    lengths = [e - s + 1 for s, e, _ in windows]
    norm = Normalize(vmin=min(lengths), vmax=max(lengths))

    for (start, end, _), n in zip(windows[:-1], lengths[:-1]):
        x, y = _profile(sweep[(sweep["start_year"] == start) & (sweep["end_year"] == end)])
        ax.plot(x, y, color=WINDOW_CMAP(norm(n)), lw=0.6, alpha=0.85, zorder=2)
    x, y = _profile(sweep[(sweep["start_year"] == REF_START) & (sweep["end_year"] == REF_END)])
    # No halo (Hannah, 2026-10-01): it hid the long windows that hug the reference
    ax.plot(x, y, color=REF_COLOUR, lw=1.4, zorder=5,
            label=window_label(REF_START, REF_END))
    ax.set_ylabel("shoreline change rate (m/yr)")
    ax.set_title(figure_title(direction), loc="left")
    _finish_axis(ax, label_towns=True)
    ax.set_ylim(-OVERLAY_Y_HALF, OVERLAY_Y_HALF)
    ax.legend(loc="lower left")
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=WINDOW_CMAP), ax=ax,
                      pad=0.01, fraction=0.03)
    cb.set_label("window length (years of record)")

    pinned = "start" if direction == "forward" else "end"
    stem = "window_profiles_overlay_{0}_from_{1}".format(direction, wc.pinned_year(direction))
    paths = fs.save(fig, Path(out_dir) / stem, close=True)
    first = "{0}–{1}".format(*windows[0][:2])
    fs.record_caption(paths[0],
        "The CoastSat shoreline change rate (OLS, the target's estimator) "
        "at every one of the island's {n} transects, south to north, fitted "
        "over each window of a nested family with the {pin} pinned at {py}: "
        "{first} through {ref}. Windows are coloured by length, light (short) "
        "to dark (long); the {ref} reference is purple. The y-axis is zoomed "
        "to ±{half:.0f} m/yr, so a short window whose rate runs past that is "
        "cut at the edge. The r of each profile against "
        "the reference is in ../../2-r_bias_rmse/.".format(
            n=sweep["transect_id"].nunique(), pin=pinned,
            py=wc.pinned_year(direction), first=first, d=direction,
            half=OVERLAY_Y_HALF,
            ref=window_label(REF_START, REF_END)))
    return paths[0]


# One panel per window
def draw_panels(sweep, corr, windows, direction, out_dir, chosen=None, years=None):
    fs.apply_style()
    drawn = (windows if years is None else
             [w for w in windows if years[0] <= w[1] - w[0] + 1 <= years[1]])
    ncol = 4
    nrow = int(np.ceil(len(drawn) / ncol))
    # A cut-down set keeps the full figure's panel height instead of stretching to a page
    height = fs.FIG_H_MAX if years is None else min(fs.FIG_H_MAX, 1.15 * nrow + 0.9)
    fig, axes = plt.subplots(nrow, ncol, figsize=fs.figsize("double", height=height),
                             sharex=True, sharey=True, layout="constrained")
    axes = np.atleast_2d(axes)
    lengths = [e - s + 1 for s, e, _ in windows]
    norm = Normalize(vmin=min(lengths), vmax=max(lengths))
    xr, yr = _profile(sweep[(sweep["start_year"] == REF_START) & (sweep["end_year"] == REF_END)])
    r_of = corr.set_index("window")["r_vs_reference"]

    for k, ax in enumerate(axes.flat):
        if k >= len(drawn):
            ax.set_visible(False)
            # The panel above an empty slot is the bottom of its column, so it carries the x-axis
            above = axes.flat[k - ncol]
            above.tick_params(labelbottom=True)
            above.set_xlabel("GIS domain", fontsize=PANEL_LABEL_PT)
            continue
        start, end, moving = drawn[k]
        is_ref = (start, end) == (REF_START, REF_END)
        if not is_ref:
            x, y = _profile(sweep[(sweep["start_year"] == start) & (sweep["end_year"] == end)])
            ax.plot(x, y, color=WINDOW_CMAP(norm(end - start + 1)), lw=0.6, zorder=2)
        ax.plot(xr, yr, color=REF_COLOUR, lw=0.9, zorder=3)
        ax.axhline(0.0, color=fs.C["INK_MUTED"], lw=0.4, zorder=1)
        ax.set_xlim(0.5, 90.5)
        ax.set_ylim(-PANEL_Y_HALF, PANEL_Y_HALF)
        ax.grid(True, axis="y", alpha=0.5)
        r = r_of.get("{0}_{1}".format(start, end), np.nan)
        title = ("{0}  ({1} yr)  reference".format(window_label(start, end), end - start + 1)
                 if is_ref else "{0}  ({1} yr)  r = {2:.2f}".format(
                     window_label(start, end), end - start + 1, r))
        ax.set_title(title, fontsize=PANEL_TITLE_PT, loc="left", color=fs.INK)
        ax.set_yticks([-5, -2.5, 0, 2.5, 5])
        ax.tick_params(labelsize=PANEL_TICK_PT)
        if chosen == (start, end):
            for spine in ax.spines.values():
                spine.set_edgecolor(CHOSEN_COLOUR)
                spine.set_linewidth(CHOSEN_LW)
                spine.set_zorder(10)
    for ax in axes[-1]:
        ax.set_xlabel("GIS domain", fontsize=PANEL_LABEL_PT)
    for ax in axes[:, 0]:
        ax.set_ylabel("m/yr", fontsize=PANEL_LABEL_PT)
    fig.suptitle(figure_title(direction), x=0.01, ha="left")

    stem = "window_profiles_panels_{0}_from_{1}".format(direction, wc.pinned_year(direction))
    if years:
        stem += "_v3"
    elif chosen:
        stem += "_v2"
    paths = fs.save(fig, Path(out_dir) / stem, close=True)
    picked = (" The chosen window, {0}, is outlined in yellow.".format(window_label(*chosen))
              if chosen else "")
    if years:
        picked += (" Only the windows {0} to {1} years long are drawn ({2} to {3}); the "
                   "full family is in the figure without _v3.".format(
                       years[0], years[1], window_label(*drawn[0][:2]),
                       window_label(*drawn[-1][:2])))
    fs.record_caption(paths[0],
        "One panel per window: the CoastSat shoreline change rate at all {n} "
        "transects (blue, darker for longer windows) against the {ref} rate "
        "(purple), on one shared y-axis zoomed to ±{half:.0f} m/yr: a short "
        "window whose rate runs past that is cut at the panel edge. Domains 1 "
        "(Cape Point) to 90 (Pea Island). Each panel names its window, its "
        "length in years of record and the alongshore Pearson r of its "
        "profile against the reference"
        + ("; the last panel is the reference alone, for comparison" if years is None
           else "") + ". The windows are nested, so r rises to 1 at the reference by "
        "construction.{picked}".format(
            n=sweep["transect_id"].nunique(), ref=window_label(REF_START, REF_END),
            half=PANEL_Y_HALF, picked=picked))
    return paths[0]


# Readme

README = """# 1-rate_profiles/{folder} — what does each window's alongshore profile look like?

Every window of the nested family ({first} … {ref}, the {pin} pinned at
{py}), drawn as a shoreline change rate profile along all {n} CoastSat
transects over the {ref} reference. Windows from **two** years, the minimum
Hannah asked for. Built {today}.

```
window_profiles_overlay_{direction}_from_{py}.png  every window over the reference (±10 m/yr)
window_profiles_panels_{direction}_from_{py}.png   one panel per window, the reference last (±5 m/yr)
{v2_line}{v3_line}window_profiles_transects.csv                     a row per transect per window (lrr, unc, n_obs, x_domain, diff)
```

How close each profile is to the reference, as r, bias and RMSE with 95%
intervals: `../../2-r_bias_rmse/`. How many years each place needs:
`../../3-settling_window/{folder}/`.

Producer:
`scripts/input_prep/5-scr/3-rates/coastsat/window_convergence/coastsat_window_profiles.py`
`--direction {direction}`.
"""


# The folder README
def write_readme(out_dir, direction, corr, n_transects):
    text = README.format(
        folder="{0}_from_{1}".format(direction, wc.pinned_year(direction)),
        ref=window_label(REF_START, REF_END),
        first=corr["window"].iloc[0].replace("_", "–"),
        pin="start" if direction == "forward" else "end",
        py=wc.pinned_year(direction), n=n_transects, direction=direction,
        today=datetime.date.today().isoformat(),
        v2_line=("window_profiles_panels_{0}_from_{1}_v2.png  the same, {2} outlined in yellow as the chosen window\n"
                 .format(direction, wc.pinned_year(direction), window_label(*CHOSEN_WINDOW[direction]))
                 if direction in CHOSEN_WINDOW else ""),
        v3_line=("window_profiles_panels_{0}_from_{1}_v3.png  the same, only the {2}-{3} yr windows\n"
                 .format(direction, wc.pinned_year(direction), *V3_YEARS[direction])
                 if direction in V3_YEARS else ""))
    (Path(out_dir) / "README.md").write_text(text, encoding="utf-8")


# One direction: profiles, figures, README; `redraw` reads the saved transects table instead of refitting
def one_direction(direction, redraw=False):
    out_dir = obs.window_profiles_dir(direction, wc.pinned_year(direction), REF_START, REF_END)
    out_dir.mkdir(parents=True, exist_ok=True)
    if redraw:
        windows = windows_for(direction)
        sweep = pd.read_csv(out_dir / "window_profiles_transects.csv")
    else:
        windows, sweep = run(direction)
        sweep.to_csv(out_dir / "window_profiles_transects.csv", index=False)
    corr = correlations(sweep, windows)
    print(draw_overlay(sweep, windows, direction, out_dir))
    print(draw_panels(sweep, corr, windows, direction, out_dir))
    # The chosen window only where this record's family contains it (2010-2026 is not in 1996-2024's)
    chosen = CHOSEN_WINDOW.get(direction)
    if chosen not in [(s, e) for s, e, _ in windows]:
        chosen = None
    if chosen:
        print(draw_panels(sweep, corr, windows, direction, out_dir, chosen=chosen))
    if direction in V3_YEARS:
        print(draw_panels(sweep, corr, windows, direction, out_dir,
                          chosen=chosen, years=V3_YEARS[direction]))
    write_readme(out_dir, direction, corr, sweep["transect_id"].nunique())
    print(corr[["window", "n_years", "r_vs_reference", "n_transects"]].to_string(index=False))


# Run: the chosen directions
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--direction", choices=("forward", "backward", "both"),
                    default="both")
    ap.add_argument("--redraw", action="store_true",
                    help="redraw from the saved transects table, no refit")
    args = ap.parse_args(argv)
    for d in (("forward", "backward") if args.direction == "both" else (args.direction,)):
        one_direction(d, redraw=args.redraw)


if __name__ == "__main__":
    main()
