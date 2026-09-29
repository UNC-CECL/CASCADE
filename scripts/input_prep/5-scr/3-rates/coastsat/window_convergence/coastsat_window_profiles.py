"""
coastsat_window_profiles.py
==============================================================================
AT WHAT WINDOW DOES THE ALONGSHORE RATE PROFILE START TO LOOK LIKE 1996-2024?

A companion to `coastsat_window_convergence.py`, drawn the other way round.
That script asks, point by point, when each transect's rate settles into a
tolerance of the reference. This one draws the WHOLE PROFILE -- rate against
alongshore position -- once per window, over the 1996-2024 reference, and
scores each window by one number: the alongshore Pearson r against the
reference, i.e. whether the pattern of erosion and accretion hotspots is right
regardless of offset (Hannah, 2026-09-29, by interview).

THE SAME NESTED FAMILIES, FROM TWO YEARS (Hannah, 2026-09-29):

    forward_from_1996/alongshore_profiles/    1996-1997, 1996-1998 ... 1996-2024
    backward_from_2024/alongshore_profiles/   2023-2024, 2022-2024 ... 1996-2024

The convergence sweep stops at five years because below that an OLS is
describing a storm cycle, not a trend. Hannah asked to see those windows
anyway, on a full y-axis: seeing how far off the short windows are is part of
the question. The minimum is therefore a SEPARATE constant here, and this
script never writes into the convergence sweep's folders or tables, so its
five-year scoring is untouched.

THE UNIT is the transect, all 906 on the island, unsmoothed. The x-axis is the
GIS domain with each domain's transects spread evenly across it in alongshore
(sorted id) order, so the village bands and the domain numbers read as in
every other alongshore figure.

WHAT r CAN AND CANNOT SAY
    The windows are NESTED: each contains the one before and the reference
    contains them all, so r goes to 1 at the reference BY CONSTRUCTION. Read
    where the curve gets there and how steadily, not whether it does. r is
    blind to offset and scale -- a window that has every hotspot in the right
    place at twice the rate scores 1.0 -- and the convergence sweep's bias and
    tolerance tables are the magnitude side of the same question.

Each window is fitted by `sweep_one` from the convergence script, which calls
the target's own `coastsat_lrr.compute_lrr`, so the estimator is the target's.
A window with fewer than MIN_OBS positions on a transect gets no fit there,
and r for that window is over the transects that have one (n in the table).

Outputs  (under obs.window_convergence_dir(direction, anchor)/alongshore_profiles/)
    window_profiles_<direction>_from_<year>.png         (a) every window over
                                                        the reference, (b) r
    window_profiles_panels_<direction>_from_<year>.png  one panel per window
    window_profiles_transects.csv                       a row per transect per window
    window_profiles_correlation.csv                     a row per window: r, n
    README.md, supporting/ (PDFs, CAPTIONS.md)

Usage
-----
    python .../coastsat_window_profiles.py                     both directions
    python .../coastsat_window_profiles.py --direction forward
==============================================================================
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

# The sibling sweep, for its loader and its per-transect fit, so the two
# products cannot fit a window differently.
sys.path.insert(0, str(Path(__file__).resolve().parent))
import coastsat_window_convergence as wc              # noqa: E402

import matplotlib                                     # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                       # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, Normalize  # noqa: E402
from matplotlib.cm import ScalarMappable              # noqa: E402

# ============================================================
# CONFIG
# ============================================================

REF_START, REF_END = wc.REF_START, wc.REF_END

# TWO YEARS, not the sweep's five (Hannah, 2026-09-29): the short windows are
# drawn so their distance from the reference can be seen. Two is the least an
# OLS through a year boundary can mean anything at all.
MIN_WINDOW_YEARS = 2

SUBFOLDER = "alongshore_profiles"

# Windows light (short) to dark (long) on the shoreline-blue ramp, coloured by
# LENGTH so both directions read the same way. The reference is the purple
# ACCENT, off the ramp, so it cannot be mistaken for a long window.
WINDOW_CMAP = LinearSegmentedColormap.from_list(
    "window_length", ["#c6dbef", "#6baed6", fs.C_1997, "#08306b"])
REF_COLOUR = fs.C["ACCENT"]


def windows_for(direction):
    """The nested family, SHORTEST FIRST, so the reference is always last."""
    if direction == "forward":
        return [(REF_START, y, y)
                for y in range(REF_START + MIN_WINDOW_YEARS - 1, REF_END + 1)]
    if direction == "backward":
        return [(y, REF_END, y)
                for y in range(REF_END - MIN_WINDOW_YEARS + 1, REF_START - 1, -1)]
    raise ValueError("direction is 'forward' or 'backward', not {0!r}".format(direction))


def window_label(start, end):
    return "{0}–{1}".format(start, end)


# ============================================================
# THE SWEEP
# ============================================================

def alongshore_x(picks):
    """Each transect's x: its domain, with the domain's transects spread
    evenly across [d - 0.5, d + 0.5] in alongshore order."""
    frame = pd.DataFrame(picks, columns=["domain_number", "transect_id"])
    rank = frame.groupby("domain_number").cumcount()
    n = frame.groupby("domain_number")["transect_id"].transform("size")
    frame["x_domain"] = frame["domain_number"] - 0.5 + (rank + 0.5) / n
    return frame


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
    sweep = sweep.merge(alongshore_x(picks)[["transect_id", "x_domain"]],
                        on="transect_id", how="left")
    ref = (sweep[(sweep["start_year"] == REF_START) & (sweep["end_year"] == REF_END)]
           .set_index("transect_id")["lrr_m_yr"])
    sweep["ref_lrr_m_yr"] = sweep["transect_id"].map(ref)
    sweep["diff_m_yr"] = sweep["lrr_m_yr"] - sweep["ref_lrr_m_yr"]
    return windows, sweep


def correlations(sweep, windows):
    """One row per window: alongshore Pearson r against the reference."""
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


# ============================================================
# FIGURES
# ============================================================

def _profile(frame):
    """x and rate for one window, in alongshore order, NaN where no fit, so
    the line breaks instead of bridging a gap."""
    frame = frame.sort_values("x_domain")
    return frame["x_domain"].to_numpy(), frame["lrr_m_yr"].to_numpy()


def _finish_axis(ax, label_towns):
    ax.axhline(0.0, color=fs.C["INK_MUTED"], lw=0.5, zorder=1)
    ax.set_xlim(0.5, 90.5)
    ax.grid(True, axis="y", alpha=0.6)
    fs.town_bands(ax, label=label_towns)


def draw_overlay(sweep, corr, windows, direction, out_dir):
    fs.apply_style()
    fig, (ax, axr) = plt.subplots(
        2, 1, figsize=fs.figsize("double", height=6.4), layout="constrained",
        gridspec_kw=dict(height_ratios=[2.2, 1.0]))
    lengths = [e - s + 1 for s, e, _ in windows]
    norm = Normalize(vmin=min(lengths), vmax=max(lengths))

    for (start, end, _), n in zip(windows[:-1], lengths[:-1]):
        x, y = _profile(sweep[(sweep["start_year"] == start) & (sweep["end_year"] == end)])
        ax.plot(x, y, color=WINDOW_CMAP(norm(n)), lw=0.6, alpha=0.85, zorder=2)
    x, y = _profile(sweep[(sweep["start_year"] == REF_START) & (sweep["end_year"] == REF_END)])
    ax.plot(x, y, color=REF_COLOUR, lw=1.4, zorder=5,
            path_effects=fs._halo(2.6), label=window_label(REF_START, REF_END))
    ax.set_ylabel("shoreline change rate (m/yr)")
    fs._title(ax, 0, "Every window over the {0} rate".format(
        window_label(REF_START, REF_END)))
    _finish_axis(ax, label_towns=True)
    ax.legend(loc="lower left")
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=WINDOW_CMAP), ax=ax,
                      pad=0.01, fraction=0.03)
    cb.set_label("window length (years of record)")

    pinned = "start" if direction == "forward" else "end"
    axr.plot(corr["n_years"], corr["r_vs_reference"], color=fs.C["LATE"], lw=1.3,
             marker="o", ms=2.8, zorder=3)
    marked = corr[corr["moving_year"] == wc.MARKED_YEAR]
    if len(marked):
        m = marked.iloc[0]
        axr.plot(m["n_years"], m["r_vs_reference"], "o", ms=6, mfc="none",
                 mec=fs.C["ADDED"], mew=1.4, zorder=4)
        axr.annotate("{0}\nr = {1:.2f}".format(
                         window_label(int(m["start_year"]), int(m["end_year"])),
                         m["r_vs_reference"]),
                     (m["n_years"], m["r_vs_reference"]), xytext=(8, -4),
                     textcoords="offset points", fontsize=7, va="top",
                     color=fs.C["ADDED"])
    axr.set_xlim(min(lengths) - 0.5, max(lengths) + 0.5)
    axr.set_ylim(min(-0.05, np.nanmin(corr["r_vs_reference"]) - 0.05), 1.05)
    axr.axhline(1.0, color=fs.C["INK_MUTED"], lw=0.5, ls=(0, (1, 2)))
    axr.grid(True, alpha=0.6)
    axr.set_xlabel("window length (years of record, {0} pinned at {1})".format(
        pinned, wc.pinned_year(direction)))
    axr.set_ylabel("alongshore r")
    fs._title(axr, 1, "Correlation of each window's profile with the reference")

    stem = "window_profiles_{0}_from_{1}".format(direction, wc.pinned_year(direction))
    paths = fs.save(fig, Path(out_dir) / stem, close=True)
    first = corr["window"].iloc[0].replace("_", "–")
    fs.record_caption(paths[0],
        "(a) The CoastSat shoreline change rate (OLS, the target's estimator) "
        "at every one of the island's {n} transects, south to north, fitted "
        "over each window of a nested family with the {pin} pinned at {py}: "
        "{first} through {ref}. Windows are coloured by length, light (short) "
        "to dark (long); the {ref} reference is purple. The y-axis holds every "
        "window, so the shortest set its range. (b) Pearson r between each "
        "window's alongshore profile and the reference, over the transects "
        "with a fit; the {mk} window is ringed. The windows are nested, so r "
        "reaches 1 at the reference by construction: read where it gets there, "
        "not whether. r ignores offset and scale.".format(
            n=sweep["transect_id"].nunique(), pin=pinned,
            py=wc.pinned_year(direction), first=first,
            ref=window_label(REF_START, REF_END),
            mk=wc.window_label(direction, wc.MARKED_YEAR)))
    return paths[0]


def draw_panels(sweep, corr, windows, direction, out_dir):
    fs.apply_style()
    drawn = windows[:-1]
    ncol = 4
    nrow = int(np.ceil(len(drawn) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=fs.figsize("double", height=fs.FIG_H_MAX),
                             sharex=True, sharey=True, layout="constrained")
    axes = np.atleast_2d(axes)
    lengths = [e - s + 1 for s, e, _ in windows]
    norm = Normalize(vmin=min(lengths), vmax=max(lengths))
    xr, yr = _profile(sweep[(sweep["start_year"] == REF_START) & (sweep["end_year"] == REF_END)])
    r_of = corr.set_index("window")["r_vs_reference"]

    for k, ax in enumerate(axes.flat):
        if k >= len(drawn):
            ax.set_visible(False)
            continue
        start, end, moving = drawn[k]
        x, y = _profile(sweep[(sweep["start_year"] == start) & (sweep["end_year"] == end)])
        ax.plot(x, y, color=WINDOW_CMAP(norm(end - start + 1)), lw=0.6, zorder=2)
        ax.plot(xr, yr, color=REF_COLOUR, lw=0.9, zorder=3)
        ax.axhline(0.0, color=fs.C["INK_MUTED"], lw=0.4, zorder=1)
        ax.set_xlim(0.5, 90.5)
        ax.grid(True, axis="y", alpha=0.5)
        r = r_of.get("{0}_{1}".format(start, end), np.nan)
        is_marked = moving == wc.MARKED_YEAR
        ax.set_title("{0}  ({1} yr)  r = {2:.2f}".format(
                         window_label(start, end), end - start + 1, r),
                     fontsize=7, loc="left",
                     color=fs.C["ADDED"] if is_marked else fs.INK,
                     fontweight="bold" if is_marked else "normal")
        ax.tick_params(labelsize=6)
    for ax in axes[-1]:
        ax.set_xlabel("GIS domain", fontsize=7)
    for ax in axes[:, 0]:
        ax.set_ylabel("m/yr", fontsize=7)

    stem = "window_profiles_panels_{0}_from_{1}".format(direction, wc.pinned_year(direction))
    paths = fs.save(fig, Path(out_dir) / stem, close=True)
    fs.record_caption(paths[0],
        "One panel per window: the CoastSat shoreline change rate at all {n} "
        "transects (blue, darker for longer windows) against the {ref} rate "
        "(purple), on one shared y-axis that holds every window, domains 1 "
        "(Cape Point) to 90 (Pea Island). Each panel names its window, its "
        "length in years of record and the alongshore Pearson r of its "
        "profile against the reference; the {mk} window is titled in amber. "
        "The windows are nested, so r rises to 1 at the reference by "
        "construction.".format(
            n=sweep["transect_id"].nunique(), ref=window_label(REF_START, REF_END),
            mk=wc.window_label(direction, wc.MARKED_YEAR)))
    return paths[0]


# ============================================================
# README
# ============================================================

README = """# {folder}/alongshore_profiles — when does the profile start to look like {ref}?

Every window of the nested family ({first} … {ref}, the {pin} pinned at
{py}), drawn as a shoreline change rate profile along all {n} CoastSat
transects over the {ref} reference, and scored by the alongshore Pearson r of
each window's profile against it. Built {today} by interview (Hannah).

The companion to the convergence sweep one folder up. That sweep scores each
transect separately against tolerances, from five years. This one draws the
whole profile from **two** years, the minimum Hannah asked for, and never
writes into the sweep's tables.

```
window_profiles_{direction}_from_{py}.png         (a) every window over the reference, (b) r per window
window_profiles_panels_{direction}_from_{py}.png  one panel per window
window_profiles_transects.csv                     a row per transect per window (lrr, unc, n_obs, x_domain, diff)
window_profiles_correlation.csv                   a row per window: r_vs_reference, n_transects
```

## r against the reference

| window | years | r |
|---|---|---|
{table}

**Read r with care.** The windows are nested, so r reaches 1 at the
reference by construction. What it shows is where it gets there and how
steadily. r measures shape only: a window with every hotspot in the right
place at the wrong magnitude still scores high. Magnitude is scored by the
bias and tolerance tables in `../all_transects/`.

Producer:
`scripts/input_prep/5-scr/3-rates/coastsat/window_convergence/coastsat_window_profiles.py`
`--direction {direction}`.
"""


def write_readme(out_dir, direction, corr, n_transects):
    lines = []
    for row in corr.itertuples(index=False):
        mark = " ← model window" if row.moving_year == wc.MARKED_YEAR else ""
        lines.append("| {0}–{1}{2} | {3} | {4:.3f} |".format(
            row.start_year, row.end_year, mark, row.n_years, row.r_vs_reference))
    text = README.format(
        folder="{0}_from_{1}".format(direction, wc.pinned_year(direction)),
        ref=window_label(REF_START, REF_END),
        first=corr["window"].iloc[0].replace("_", "–"),
        pin="start" if direction == "forward" else "end",
        py=wc.pinned_year(direction), n=n_transects, direction=direction,
        today=datetime.date.today().isoformat(), table="\n".join(lines))
    (Path(out_dir) / "README.md").write_text(text, encoding="utf-8")


# ============================================================
# MAIN
# ============================================================

def one_direction(direction):
    out_dir = obs.window_convergence_dir(direction, wc.pinned_year(direction),
                                         ref_start=REF_START, ref_end=REF_END) / SUBFOLDER
    out_dir.mkdir(parents=True, exist_ok=True)
    windows, sweep = run(direction)
    corr = correlations(sweep, windows)
    sweep.to_csv(out_dir / "window_profiles_transects.csv", index=False)
    corr.to_csv(out_dir / "window_profiles_correlation.csv", index=False)
    print(draw_overlay(sweep, corr, windows, direction, out_dir))
    print(draw_panels(sweep, corr, windows, direction, out_dir))
    write_readme(out_dir, direction, corr, sweep["transect_id"].nunique())
    print(corr[["window", "n_years", "r_vs_reference", "n_transects"]].to_string(index=False))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[3])
    ap.add_argument("--direction", choices=("forward", "backward", "both"),
                    default="both")
    args = ap.parse_args(argv)
    for d in (("forward", "backward") if args.direction == "both" else (args.direction,)):
        one_direction(d)


if __name__ == "__main__":
    main()
