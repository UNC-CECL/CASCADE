"""
Two SOURCES for one start year: the dune line against the shoreline
==============================================================================
Compares the island offset built from the digitised DUNE line with the one
built from the CoastSat SHORELINE, for one hindcast start. Written 2026-09-22,
the day the shoreline source was added.

NOT THE SAME QUESTION AS HAT_compare_offset_versions.py
    That script compares two VERSIONS of one source: the same feature
    re-digitised, where the difference is a correction and the interesting
    number is how many domains moved. This one compares two SOURCES: two
    DIFFERENT FEATURES on the island, where the difference is the beach
    between them and is supposed to be there. Different question, different
    framing, different caption -- so a sibling script rather than a flag.

THE TRAP THIS FIGURE EXISTS TO SHOW
    Correlate the two profiles and you get r = 1.0000, which looks like the
    two features agreeing about the island. They do not. Both are dominated by
    the same ~6.2 km of cape curvature, against which the 17 m standard
    deviation of their difference is 0.3%. So the figure states the
    correlation nowhere and shows the gap instead.

WHICH FRAME IS DRAWN, AND WHY IT HAD TO CHANGE
    Both profiles are drawn as the STATION FROM THE SHARED OFFSHORE DATUM, not
    as the min-zeroed offset the model is handed.

    The figure was drawn in the model frame until 2026-09-22, when Hannah
    asked the question that breaks it: "shouldn't the dune always be behind
    the shoreline?" It should, and on the ground it is -- the shoreline is
    seaward of the dune line in 90 of 90 domains, by 17-97 m. But the two
    model files are each zeroed on their OWN most seaward domain, and those
    minima are 45.6 m apart (dune 1953.198, shoreline 1907.607), so

        model_diff = -beach_width + 45.6 m

    and the dune line comes out apparently seaward wherever the beach is
    narrower than 45.6 m. That was 69 of 90 domains -- exactly the 69 the old
    figure shaded "dune line seaward", a physically impossible claim generated
    entirely by the zeroing.

    The datum frame has no such constant, and it is the SAME SHAPE (a model
    offset is this minus the build's own minimum). So the profiles are
    unchanged, the band between them is the beach, and it is on the correct
    side of the island everywhere. What the model reads is still one
    subtraction away, and the CSV holds both frames.

HOW IT IS DRAWN, AND WHAT THAT COSTS
    Six consecutive sections on a 2 x 3 grid, each a VERTICAL strip:
    alongshore up the page, cross-shore across it, south at the bottom, the
    axis inverted so the ocean is on the right. The panels read as the island.
    Why a GRID and not a row of six: see SECTIONS below -- columns take width
    from the very axis the gap is measured on.

    These are the ABSOLUTE profiles, and the beach is a small fraction of the
    cross-shore range they span, so the band is thin. Removing a smooth trend
    from both sources would open it up (14-43% of a panel rather than 2-3%),
    and was built and then taken back out (Hannah, 2026-09-22: "I actually
    don't like the detrend, I want to see the original shoreline shape"). The
    shape is the point; each panel therefore states its own beach range as a
    number.

OUTPUT   2-brie-offset/<year>/comparisons/<a>_vs_<b>/
    offset_<year>_<a>_vs_<b>.csv        per domain, both frames, columns named
                                        for the SOURCE not for a version
    offset_<year>_<a>_vs_<b>.png/.pdf   the two profiles as vertical strips of
                                        island, caption in CAPTIONS.md
    README.md                           what the folder is

    A comparison is neither a version nor a source, so it lands in
    comparisons/ rather than inside either build
    (hat_topo_version.offset_comparison_dir).

USAGE
    python compare_offset_sources.py --year 1996
    python compare_offset_sources.py --year 1996 --a duneline --b shoreline
==============================================================================
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

PROJECT_ROOT = next(_p for _p in Path(__file__).resolve().parents
                    if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from site_layer import hat_topo_version as _tv  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C, FIG_W_DOUBLE, INK_MUTED, _title, apply_style, caption, save)

# HOW EACH SOURCE IS NAMED AND DRAWN.
#
# The house RdBu pair (C_1984 red / C_1997 blue) means EARLIER and LATER
# vintage, and these two are not a time order -- drawing them red and blue
# would say the dune line came first. Blue is already the CoastSat colour in
# this style sheet, so the shoreline keeps it and the dune line takes the
# accent purple.
SOURCE_STYLE = {
    "duneline": {
        "label": "dune line",
        "colour": C["ACCENT"],
        "fill": C["ACCENT_FILL"],
        "feature": "dune line digitised from aerial imagery",
    },
    "shoreline": {
        "label": "CoastSat shoreline",
        "colour": C["LATE"],
        "fill": C["LATE_FILL"],
        "feature": "mean satellite shoreline over an averaging window",
    },
}


# The island in sixths, drawn as vertical strips on a 2 x 3 GRID (2026-09-22).
#
# THE GRID IS THE POINT, not the section count. These panels are columns, so
# every extra one alongside takes width away from the offset axis -- the axis
# the gap is measured on -- faster than the tighter zoom gives back. Measured,
# as the widest gap actually rendered on paper at 300 dpi:
#
#     sections, across   panel width   gap on paper
#        3                  2.07 in       13 px
#        4                  1.52 in       10 px
#        6                  0.96 in        7 px   <- MORE panels, LESS gap
#       10                  0.52 in        6 px
#
# Splitting into ROWS instead gives the width back, so the zoom is kept and
# the offset axis is not squeezed:
#
#     6 as 2 x 3          2.07 in       15 px
#     9 as 3 x 3          2.07 in       23 px
#
# 6 on a 2 x 3 grid is the chosen point: half again the separation of 4 across,
# 15 domains a panel, and panels still tall enough (3.9 in) to read as strips
# of coast.
SECTIONS = ((1, 15), (16, 30), (31, 45), (46, 60), (61, 75), (76, 90))
GRID_COLS = 3


def _town_bands_alongshore_y(ax, lo, hi):
    """Village spans, for a panel whose ALONGSHORE axis is the vertical one.

    hat_figure_style.town_bands draws the same thing against a horizontal
    alongshore axis and has no orientation switch, so this is its axhspan
    twin. The spans themselves still come from the one owner,
    hatteras_site_config.HATTERAS_ANNOTATIONS -- only the axis differs.
    """
    try:
        from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS
        spans = HATTERAS_ANNOTATIONS.town_spans
    except ImportError:
        return
    for name, (s_lo, s_hi) in spans.items():
        if s_hi + 0.5 < lo - 0.5 or s_lo - 0.5 > hi + 0.5:
            continue
        ax.axhspan(s_lo - 0.5, s_hi + 0.5, color="0.94", lw=0, zorder=0)
        mid = (max(s_lo - 0.5, lo - 0.5) + min(s_hi + 0.5, hi + 0.5)) / 2
        ax.text(0.015, mid, name, transform=ax.get_yaxis_transform(),
                ha="left", va="center", fontsize=6.5, color=INK_MUTED,
                zorder=1, clip_on=True)


def _window_gap_months(year):
    """Months between the centre of the shoreline averaging window and the
    dune line's survey date, or None if either is unavailable.

    Computed, never typed: it was typed once, as "about nine months", and it
    is 15.4 (found 2026-09-22 when the figure asked for it). The survey date
    has one owner, duneline_endpoint.survey_date, so this asks it.
    """
    try:
        import datetime as dt
        sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "input_prep" / "5-scr" / "lib"))
        import scr_paths  # noqa: F401  (puts the 5-scr modules on sys.path)
        from duneline_endpoint import survey_date
        start, end = _tv.shoreline_window_for_year(year)
        centre = (dt.date(start, 1, 1)
                  + (dt.date(end, 12, 31) - dt.date(start, 1, 1)) / 2)
        surveyed, _assumed = survey_date(_tv.dune_line_for_year(year))
        return abs((surveyed - centre).days) / 30.44
    except Exception:
        return None


def _gap_clause(year):
    """The time between the two observations, spelled for a caption."""
    months = _window_gap_months(year)
    return ("" if months is None else
            ", plus the {0:.0f} months between the centre of the shoreline "
            "window and the dune line's survey".format(months))


def _vintage_label(source, year):
    """What the source IS for this start year, spelled for a legend: a dune
    line carries an imagery vintage, a shoreline carries a window."""
    if source == "duneline":
        return "dune line ({0} imagery)".format(_tv.dune_line_for_year(year))
    window = _tv.shoreline_window_for_year(year)
    return "CoastSat shoreline ({0}–{1} mean)".format(*window)


def _unpadded(year, source):
    """The 90-domain file the model would read from one source's CURRENT
    build, zeroed on that build's own most seaward domain."""
    path = _tv.offset_file(year, "unpadded", source=source)
    if not path.is_file():
        sys.exit("no {0} build for {1}: {2} is missing".format(source, year, path))
    return pd.read_csv(path).set_index("Domain_ID")[str(year)]


def _raw_domain_means(path):
    """Per-domain mean station from the shared offshore datum. One row per
    transect first, so a domain with more transects does not weight twice."""
    raw = pd.read_csv(path)
    per_transect = raw.drop_duplicates(["domain_id", "LineID"])
    return per_transect.groupby("domain_id")["ORIG_LEN"].mean()


def _raw_file(year, source):
    return (_tv.dune_raw_file_for_year(year) if source == "duneline"
            else _tv.shoreline_raw_file_for_year(year))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--year", type=int, default=1996)
    ap.add_argument("--a", default="duneline", choices=tuple(SOURCE_STYLE),
                    help="the reference source (default duneline)")
    ap.add_argument("--b", default="shoreline", choices=tuple(SOURCE_STYLE),
                    help="the source compared against it (default shoreline)")
    args = ap.parse_args(argv)
    if args.a == args.b:
        ap.error("--a and --b name the same source")
    year, a, b = args.year, args.a, args.b
    lab_a, lab_b = _vintage_label(a, year), _vintage_label(b, year)

    # ---- the numbers ----------------------------------------------------- #
    ma, mb = _unpadded(year, a), _unpadded(year, b)
    out = pd.DataFrame({"model_{0}_m".format(a): ma, "model_{0}_m".format(b): mb})
    out["model_diff_m"] = mb - ma

    ra = _raw_domain_means(_raw_file(year, a))
    rb = _raw_domain_means(_raw_file(year, b))
    out["datum_{0}_m".format(a)] = ra.reindex(out.index)
    out["datum_{0}_m".format(b)] = rb.reindex(out.index)
    # Stations grow LANDWARD from the offshore datum, so a - b is positive
    # where b lies seaward of a. With a=duneline and b=shoreline that is the
    # beach width.
    out["seaward_gap_m"] = out["datum_{0}_m".format(a)] - out["datum_{0}_m".format(b)]
    out.index.name = "gis_domain"

    gap, mdiff = out["seaward_gap_m"], out["model_diff_m"]
    baseline_shift = float(ra.min() - rb.min())
    print("{0}: {1} vs {2}".format(year, a, b))
    print("  fixed datum -- {0} seaward of {1} by: mean {2:+.1f} m, median {3:+.1f}, "
          "sd {4:.1f}, range {5:+.1f} .. {6:+.1f}".format(
              b, a, gap.mean(), gap.median(), gap.std(), gap.min(), gap.max()))
    print("  {0} of {1} domains have {2} seaward".format(
        int((gap > 0).sum()), len(gap), b))
    print("  baseline gap ({0} min - {1} min): {2:+.1f} m".format(a, b, baseline_shift))
    print("  model frame -- {0} minus {1}: mean {2:+.1f} m, sd {3:.1f}, "
          "range {4:+.1f} .. {5:+.1f}".format(
              b, a, mdiff.mean(), mdiff.std(), mdiff.min(), mdiff.max()))
    print("  (model frame = -(seaward gap) + {0:+.1f} m, which is why the sign "
          "flips)".format(baseline_shift))

    out_dir = _tv.offset_comparison_dir(year, "{0}_vs_{1}".format(a, b))
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = "offset_{0}_{1}_vs_{2}".format(year, a, b)
    out.to_csv(out_dir / "{0}.csv".format(stem), float_format="%.2f")

    # ---- figure ---------------------------------------------------------- #
    # FOUR SECTIONS, DRAWN VERTICALLY, ABSOLUTE OFFSETS (2026-09-22, Hannah:
    # "I want to see the original shoreline shape ... should you do vertical
    # instead like the actual island shape").
    #
    # Alongshore runs UP the page and the offset across it, so the panels read
    # as four consecutive strips of Hatteras with south at the bottom. The x
    # axis is INVERTED: an offset grows landward (the 100 m transects run due
    # west from the offshore datum), so inverting it puts the ocean on the
    # right, where it is.
    #
    # No detrending: these are the profiles as the model reads them. That
    # costs legibility and it is worth being honest about how much -- a
    # quarter of the island still spans 800-2200 m of offset, so the widest
    # gap between the two sources is 1.3-3.3% of a panel's width. The filled
    # band is what carries it; in the flatter sections it is a sliver.
    col_a, col_b = SOURCE_STYLE[a]["colour"], SOURCE_STYLE[b]["colour"]
    fill_b = SOURCE_STYLE[b]["fill"]
    # DRAWN IN THE FIXED-DATUM FRAME, not the model frame (2026-09-22, Hannah:
    # "shouldn't the dune always be behind the shoreline?").
    #
    # It should, and on the ground it is: the shoreline is seaward of the dune
    # line in 90 of 90 domains. But the MODEL files are each zeroed on their
    # own most seaward domain, and those two minima are 45.6 m apart, so
    # differencing them puts the dune line apparently seaward wherever the
    # beach is narrower than 45.6 m -- which is 69 of 90 domains. The earlier
    # version of this figure drew exactly that and labelled it "dune line
    # seaward", a physically impossible claim produced entirely by the zeroing.
    #
    # The station from the shared offshore datum carries no such constant, and
    # it is the SAME SHAPE: a model offset is this minus the build's own
    # minimum. So nothing about the profiles is lost, the band between them is
    # the beach, and it is on the correct side everywhere.
    ya, yb = out["datum_{0}_m".format(a)], out["datum_{0}_m".format(b)]

    apply_style()
    n_rows = int(np.ceil(len(SECTIONS) / GRID_COLS))
    fig, axes = plt.subplots(n_rows, GRID_COLS,
                             figsize=(FIG_W_DOUBLE, 3.9 * n_rows + 1.1),
                             constrained_layout=True)
    axes = np.atleast_1d(axes).ravel()
    for spare in axes[len(SECTIONS):]:
        spare.set_visible(False)

    for i, (lo, hi) in enumerate(SECTIONS):
        ax = axes[i]
        sl = (out.index >= lo) & (out.index <= hi)
        dom = out.index.to_numpy()[sl]
        va, vb = ya.to_numpy()[sl], yb.to_numpy()[sl]

        _town_bands_alongshore_y(ax, lo, hi)

        # One band, one colour, one direction: the station grows LANDWARD from
        # the datum and the shoreline's is always the smaller, so the band is
        # always the beach. It needs no second colour for a sign that cannot
        # change.
        ax.fill_betweenx(dom, va, vb, color=fill_b, lw=0, zorder=2,
                         label="beach (shoreline to dune line)")
        ax.plot(va, dom, color=col_a, lw=1.3, zorder=4, label=lab_a)
        ax.plot(vb, dom, color=col_b, lw=1.3, ls=(0, (4, 2)), zorder=5,
                label=lab_b)

        ax.set_ylim(lo - 0.5, hi + 0.5)
        ax.invert_xaxis()          # landward left, ocean right
        ax.xaxis.set_major_locator(mticker.MaxNLocator(3))
        # A domain is a count, so its ticks are whole numbers; the default
        # locator offered 42.5 and 45.0.
        ax.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))
        ax.tick_params(axis="x", labelsize=7)

        beach = va - vb          # + = shoreline seaward of the dune line
        k = int(np.argmax(beach))
        _title(ax, i, "GIS {0}–{1}".format(lo, hi))
        # The beach as a number, inside the panel: at this scale the eye
        # cannot measure the band, and the title has no room for it.
        ax.annotate("beach {0:.0f}–{1:.0f} m\nwidest at GIS {2}".format(
                        beach.min(), beach.max(), dom[k]),
                    xy=(0.5, 0.008), xycoords="axes fraction", ha="center",
                    va="bottom", fontsize=6.5, color=INK_MUTED,
                    bbox=dict(facecolor="white", alpha=0.85, edgecolor="none",
                              boxstyle="square,pad=0.2"))

    fig.supylabel("GIS domain (south → north)", fontsize=9)
    fig.supxlabel("Distance from the offshore datum (m)   —   "
                  "landward ←   |   → ocean", fontsize=9)
    # ABOVE the panels, not below: at the bottom the legend and the shared
    # offset label are both "outside lower centre" and constrained_layout
    # stacks them on top of each other.
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside upper center", ncol=3, fontsize=7,
               frameon=False)

    caption(fig, (
        "The two features the {year} island offset can be built from — the {fa}, "
        "and the {fb} — in {nsec} consecutive sections of Hatteras on a "
        "{nrow} x {ncol} grid. Alongshore runs up the page and cross-shore across "
        "it, so each panel is a strip of the island with south at the bottom and "
        "the ocean on the right. Both are drawn as the distance from the SHARED "
        "offshore datum along the 100 m transects, averaged per 500 m domain, so "
        "the two are in one frame and the band between them is the beach: the "
        "shoreline lies seaward of the dune line in all {n} domains, by {gm:+.1f} m "
        "on average and {glo:.0f}–{ghi:.0f} m across the island{gapclause}. Each "
        "panel gives its own beach range, because at this scale the eye cannot "
        "measure the band. "
        "The profiles are the same SHAPE the model reads, but not the same "
        "numbers: each build is handed to the model zeroed on its own most seaward "
        "domain, and those two minima are {shift:.1f} m apart. Differencing the two "
        "model files therefore subtracts that constant and reports the dune line "
        "as the seaward feature wherever the beach is narrower than it — which is "
        "{nflip} of {n} domains, and is an artefact of the zeroing, not a thing "
        "that happens on this island. The per-domain numbers, in both frames, are "
        "in {stem}.csv.").format(
            year=year, fa=SOURCE_STYLE[a]["feature"], fb=SOURCE_STYLE[b]["feature"],
            nsec=len(SECTIONS), nrow=n_rows, ncol=GRID_COLS,
            n=len(gap), gm=gap.mean(), glo=gap.min(), ghi=gap.max(),
            gapclause=_gap_clause(year), shift=abs(baseline_shift),
            nflip=int((mdiff > 0).sum()), stem=stem))
    paths = save(fig, out_dir / stem, close=True)
    print("  wrote {0}".format(out_dir / (stem + ".csv")))
    for p in paths:
        print("  wrote {0}".format(p))

    _write_readme(out_dir, year, a, b, lab_a, lab_b, gap, mdiff, baseline_shift, stem)
    print("  wrote {0}".format(out_dir / "README.md"))


def _write_readme(out_dir, year, a, b, lab_a, lab_b, gap, mdiff, shift, stem):
    (out_dir / "README.md").write_text("""# {year} island offset: {a} vs {b}

Two builds of the **same** {year} island offset, from two **different
features** on the island:

| source | what it is | build |
|---|---|---|
| `{a}` | {lab_a} | `../../v1/` |
| `{b}` | {lab_b} | `../../{b}/v1/` |

Written by `scripts/input_prep/2-brie-offset/2-figures/compare_offset_sources.py`.
This is a comparison, not a build: nothing here is read by a model run.

## What it found

| | |
|---|---|
| {b} seaward of {a} (fixed datum) | mean {gm:+.1f} m, median {gmed:+.1f}, sd {gsd:.1f}, range {glo:+.1f} to {ghi:+.1f} |
| domains with {b} seaward | {npos} of {n} |
| gap between the two zeroing baselines | {shift:+.1f} m |
| difference in the model frame | mean {mm:+.1f} m, sd {msd:.1f}, range {mlo:+.1f} to {mhi:+.1f} |

## Read the correlation, not the way round you expect

The two profiles correlate at **r = 1.0000**. That is not the dune line and
the shoreline agreeing about the beach — it is both of them being dominated by
the same ~6.2 km of cape curvature, against which the {msd:.0f} m sd of their
difference is 0.3%. **Difference these profiles; never correlate them.**

## The figure

Six consecutive sections of the island on a **2 × 3 grid**, each a **vertical
strip**: alongshore up the page, cross-shore across it, south at the bottom,
ocean on the right. The panels read as the island.

Both profiles are drawn as the **distance from the shared offshore datum**,
not as the min-zeroed offset the model is handed — see the next section for
why that matters. They are the **absolute** profiles, the shape each source
gives, not a residual.

## The dune line is behind the shoreline, and the figure has to say so

On the ground it always is: the shoreline is seaward of the dune line in
**90 of 90 domains**, by 17–97 m. The band between the two lines is that beach.

This figure was drawn in the **model frame** until 2026-09-22, and in that
frame it was not true. Each build is zeroed on its own most seaward domain,
and those two minima are **45.6 m** apart, so

```
model_diff = -(beach width) + 45.6 m
```

and the dune line comes out *apparently seaward* wherever the beach is
narrower than 45.6 m — **69 of 90 domains**, which the old version shaded as
"dune line seaward". That was an artefact of the zeroing, not anything that
happens on this island.

The datum frame carries no such constant and is the **same shape** (a model
offset is this minus the build's own minimum), so nothing about the profiles
was lost by moving to it. Both frames are in the CSV.

**The lesson generalises:** never difference two min-zeroed offset files and
read the sign as physical. `2-brie-offset/README.md` has warned about this
since before the shoreline source existed; this is what it looks like when it
bites.

## Why a grid, and not just more panels

These panels are columns, so each extra one **alongside** takes width from the
offset axis — the axis the gap is measured on — faster than the tighter zoom
gives back. Measured as the widest gap actually rendered, at 300 dpi:

| layout | panel width | gap on paper |
|---|---|---|
| 3 across | 2.07 in | 13 px |
| 4 across | 1.52 in | 10 px |
| **6 across** | 0.96 in | **7 px** — more panels, *less* gap |
| 10 across | 0.52 in | 6 px |

Splitting into **rows** gives the width back, so the zoom is kept and the
offset axis is not squeezed:

| layout | panel width | gap on paper |
|---|---|---|
| **6 as 2 × 3** (this figure) | 2.07 in | **15 px** |
| 9 as 3 × 3 | 2.07 in | 23 px |

Removing a smooth trend from both sources would do better still (14–43% of the
panel), and was built and then taken back out — the absolute shape is the
point. So the gap is carried by the filled band and, where the band is a
sliver, by the **widest-gap number written into each panel**. The per-domain
numbers, in both frames, are in the CSV.

## Two frames, and the sign flip between them

`ORIG_LEN` grows **landward** from the shared offshore datum, so in the
fixed-datum frame `{a} − {b}` is positive where {b} is the more seaward
feature — a beach width. The model frame is not that: each build is zeroed on
its own most seaward domain, so differencing the two subtracts the {shift:+.1f} m
gap between those baselines and flips the sign,

```
model_diff = -(seaward gap) + {shift:+.1f} m
```

which is why a mean beach width of {gm:+.1f} m appears in the model frame as
{mm:+.1f} m. **The band in the figure is the model-frame gap, so it is not a
beach width.** The `seaward_gap_m` column of the CSV is.

## Files

| file | what it is |
|---|---|
| `{stem}.csv` | per domain: both sources in both frames, columns named for the source |
| `{stem}.png` | the three-panel figure (PDF and caption under `supporting/`) |

## Rebuild

```
python scripts/input_prep/2-brie-offset/2-figures/compare_offset_sources.py --year {year}
```

Both sources resolve through their `CURRENT`, so this re-reads whatever each
source currently points at rather than the builds that were current on the day
it was written.
""".format(year=year, a=a, b=b, lab_a=lab_a, lab_b=lab_b,
           gm=gap.mean(), gmed=gap.median(), gsd=gap.std(),
           glo=gap.min(), ghi=gap.max(),
           npos=int((gap > 0).sum()), n=len(gap), shift=shift,
           mm=mdiff.mean(), msd=mdiff.std(), mlo=mdiff.min(), mhi=mdiff.max(),
           stem=stem), encoding="utf-8")


if __name__ == "__main__":
    main()
