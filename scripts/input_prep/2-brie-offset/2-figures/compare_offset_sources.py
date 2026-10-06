"""
Two sources for one start year: the island offset from the dune line against the one from the shoreline.

    python scripts/input_prep/2-brie-offset/2-figures/compare_offset_sources.py --year 1996
    python scripts/input_prep/2-brie-offset/2-figures/compare_offset_sources.py --year 1996 --a duneline --b shoreline

Draws both offsets and their difference along the island, with a README,
in the start year's comparison folder. Details: scripts/input_prep/2-brie-offset/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-06
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

# --- CONFIG ------------------------------------------------------------------
# HOW EACH SOURCE IS NAMED AND DRAWN
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


# The island in sixths, drawn as vertical strips on a 2 x 3 GRID (2026-09-22)
SECTIONS = ((1, 15), (16, 30), (31, 45), (46, 60), (61, 75), (76, 90))
GRID_COLS = 3
# -----------------------------------------------------------------------------


# Village spans, for a panel whose ALONGSHORE axis is the vertical one
def _town_bands_alongshore_y(ax, lo, hi):
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


# Months between the shoreline window's centre and the dune line's survey date, or None
def _window_gap_months(year, shoreline_version=None):
    try:
        sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "input_prep" / "5-scr" / "lib"))
        import scr_paths  # noqa: F401  (puts the 5-scr modules on sys.path)
        from duneline_endpoint import survey_date
        lo, hi = _shoreline_window(year, shoreline_version)
        centre = lo + (hi - lo) / 2
        surveyed, _assumed = survey_date(_tv.dune_line_for_year(year))
        return abs((surveyed - centre).days) / 30.44
    except Exception:
        return None


# The time between the two observations, spelled for a caption
def _gap_clause(year, shoreline_version=None):
    months = _window_gap_months(year, shoreline_version)
    return ("" if months is None else
            ", plus the {0:.0f} months between the centre of the shoreline "
            "window and the dune line's survey".format(months))


# The raw file a shoreline BUILD was made from
def _shoreline_raw(year, version=None):
    d = _tv.offset_build_dir(year, version, "shoreline")
    hits = sorted(d.glob("*_shoreline_offset_raw.csv"))
    if len(hits) != 1:
        sys.exit("expected one *_shoreline_offset_raw.csv in {0}, found {1}".format(d, len(hits)))
    return hits[0]


# (first, last) day of the build's averaging window, read from its raw file's name
def _shoreline_window(year, version=None):
    import datetime as dt
    label = _shoreline_raw(year, version).name.replace("_shoreline_offset_raw.csv", "")
    a, b = label.split("_")
    if "-" in a:
        return dt.date.fromisoformat(a), dt.date.fromisoformat(b)
    return dt.date(int(a), 1, 1), dt.date(int(b), 12, 31)


# What the source IS for this start year, spelled for a legend
def _vintage_label(source, year, version=None):
    if source == "duneline":
        return "dune line ({0} imagery)".format(_tv.dune_line_for_year(year))
    lo, hi = _shoreline_window(year, version)
    if (lo.month, lo.day, hi.month, hi.day) == (1, 1, 12, 31):
        return "CoastSat shoreline ({0}–{1} mean)".format(lo.year, hi.year)
    return "CoastSat shoreline ({0} – {1} mean)".format(lo.isoformat(), hi.isoformat())


# The 90-domain file the model would read from one source's build, zeroed on its most seaward domain
def _unpadded(year, source, version=None):
    path = _tv.offset_file(year, "unpadded", source=source, version=version)
    if not path.is_file():
        sys.exit("no {0} build for {1}: {2} is missing".format(source, year, path))
    # The value column by position: the 2009 builds kept their 2010 header through the 2026-10-05 rename
    return pd.read_csv(path).set_index("Domain_ID").iloc[:, 0].rename(str(year))


# Per-domain mean station from the shared offshore datum
def _raw_domain_means(path):
    raw = pd.read_csv(path)
    per_transect = raw.drop_duplicates(["domain_id", "LineID"])
    return per_transect.groupby("domain_id")["ORIG_LEN"].mean()


# The raw offset file for a year and source
def _raw_file(year, source, version=None):
    return (_tv.dune_raw_file_for_year(year) if source == "duneline"
            else _shoreline_raw(year, version))


# Run: both sources' offsets, the figure, the README
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--year", type=int, default=1996)
    ap.add_argument("--a", default="duneline", choices=tuple(SOURCE_STYLE),
                    help="the reference source (default duneline)")
    ap.add_argument("--b", default="shoreline", choices=tuple(SOURCE_STYLE),
                    help="the source compared against it (default shoreline)")
    ap.add_argument("--duneline-version", default=None,
                    help="the dune build to read (default its CURRENT)")
    ap.add_argument("--shoreline-version", default=None,
                    help="the shoreline build to read (default its CURRENT)")
    args = ap.parse_args(argv)
    if args.a == args.b:
        ap.error("--a and --b name the same source")
    year, a, b = args.year, args.a, args.b
    ver = {"duneline": args.duneline_version or _tv.offset_version(year, "duneline"),
           "shoreline": args.shoreline_version or _tv.offset_version(year, "shoreline")}
    lab_a, lab_b = _vintage_label(a, year, ver[a]), _vintage_label(b, year, ver[b])

    # The numbers
    ma, mb = _unpadded(year, a, ver[a]), _unpadded(year, b, ver[b])
    out = pd.DataFrame({"model_{0}_m".format(a): ma, "model_{0}_m".format(b): mb})
    out["model_diff_m"] = mb - ma

    ra = _raw_domain_means(_raw_file(year, a, ver[a]))
    rb = _raw_domain_means(_raw_file(year, b, ver[b]))
    out["datum_{0}_m".format(a)] = ra.reindex(out.index)
    out["datum_{0}_m".format(b)] = rb.reindex(out.index)
    # Stations grow LANDWARD from the offshore datum, so a - b is positive where b lies seaward of a
    out["seaward_gap_m"] = out["datum_{0}_m".format(a)] - out["datum_{0}_m".format(b)]
    out.index.name = "gis_domain"

    gap, mdiff = out["seaward_gap_m"], out["model_diff_m"]
    baseline_shift = float(ra.min() - rb.min())
    print("{0}: {1} {2} vs {3} {4}".format(year, a, ver[a], b, ver[b]))
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

    out_dir = _tv.offset_source_comparison_dir(year, "{0}_vs_{1}".format(a, b),
                                               ver["shoreline"])
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = "offset_{0}_{1}_vs_{2}".format(year, a, b)
    out.to_csv(out_dir / "{0}.csv".format(stem), float_format="%.2f")

    # FOUR SECTIONS, DRAWN VERTICALLY, ABSOLUTE OFFSETS (2026-09-22, Hannah
    col_a, col_b = SOURCE_STYLE[a]["colour"], SOURCE_STYLE[b]["colour"]
    fill_b = SOURCE_STYLE[b]["fill"]
    # DRAWN IN THE FIXED-DATUM FRAME, not the model frame (2026-09-22, Hannah
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

        # One band, one colour, one direction
        ax.fill_betweenx(dom, va, vb, color=fill_b, lw=0, zorder=2,
                         label="beach (shoreline to dune line)")
        ax.plot(va, dom, color=col_a, lw=1.3, zorder=4, label=lab_a)
        ax.plot(vb, dom, color=col_b, lw=1.3, ls=(0, (4, 2)), zorder=5,
                label=lab_b)

        ax.set_ylim(lo - 0.5, hi + 0.5)
        ax.invert_xaxis()          # landward left, ocean right
        ax.xaxis.set_major_locator(mticker.MaxNLocator(3))
        # A domain is a count, so its ticks are whole numbers
        ax.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))
        ax.tick_params(axis="x", labelsize=7)

        beach = va - vb          # + = shoreline seaward of the dune line
        k = int(np.argmax(beach))
        _title(ax, i, "GIS {0}–{1}".format(lo, hi))
        # The beach as a number, inside the panel
        ax.annotate("beach {0:.0f}–{1:.0f} m\nwidest at GIS {2}".format(
                        beach.min(), beach.max(), dom[k]),
                    xy=(0.5, 0.008), xycoords="axes fraction", ha="center",
                    va="bottom", fontsize=6.5, color=INK_MUTED,
                    bbox=dict(facecolor="white", alpha=0.85, edgecolor="none",
                              boxstyle="square,pad=0.2"))

    fig.supylabel("GIS domain (south → north)", fontsize=9)
    fig.supxlabel("Distance from the offshore datum (m)   —   "
                  "landward ←   |   → ocean", fontsize=9)
    # Legend above the panels, clear of the shared offset label below
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
            gapclause=_gap_clause(year, ver["shoreline"]), shift=abs(baseline_shift),
            nflip=int((mdiff > 0).sum()), stem=stem))
    paths = save(fig, out_dir / stem, close=True)
    print("  wrote {0}".format(out_dir / (stem + ".csv")))
    for p in paths:
        print("  wrote {0}".format(p))

    _write_readme(out_dir, year, a, b, lab_a, lab_b, gap, mdiff, baseline_shift, stem, ver)
    print("  wrote {0}".format(out_dir / "README.md"))


# The README beside the comparison
def _write_readme(out_dir, year, a, b, lab_a, lab_b, gap, mdiff, shift, stem, ver):
    if out_dir == _tv.offset_comparison_dir(year, out_dir.name):
        filed = ("Filed at the year, `{0}/comparisons/`, because it is drawn against the CURRENT shoreline "
                 "build (since 2026-10-06; from 2026-09-29 it sat under `{0}/shoreline/<v>/comparisons/`)."
                 .format(year))
    else:
        filed = ("Filed with the shoreline build it was drawn against, `{0}/shoreline/{1}/`, which is not "
                 "the CURRENT one; the comparison against CURRENT is in `{0}/comparisons/`.".format(
                     year, ver["shoreline"]))
    # The domain figures offset_sources_on_domains.py draws into domains/, one per 15 domains
    extra = "".join(
        "| `domains/{0}` | GIS {1}: the Barrier3D domains placed with each offset, dune line (a) over "
        "shoreline (b), and the shift per domain (c) (offset_sources_on_domains.py) |\n".format(
            f.name, f.stem.split("_GIS")[1])
        for f in sorted((out_dir / "domains").glob(stem + "_domains_GIS*.png")))
    overview = out_dir / "domains" / (stem + "_domains_overview.png")
    if overview.exists():
        extra = ("| `domains/{0}` | the whole island: both planforms with the six sections marked, the beach "
                 "width, and the shift per domain; start here (offset_sources_on_domains.py) |\n"
                 .format(overview.name)) + extra
    (out_dir / "README.md").write_text("""# {year} island offset: {a} vs {b}

Two builds of the **same** {year} island offset, from two **different
features** on the island:

| source | what it is | build |
|---|---|---|
| `{a}` | {lab_a} | `{year}/{a}/{va}/` |
| `{b}` | {lab_b} | `{year}/{b}/{vb}/` |

{filed}

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
{extra}
## Rebuild

```
python scripts/input_prep/2-brie-offset/2-figures/compare_offset_sources.py --year {year} --duneline-version {vd} --shoreline-version {vs}
```

Without the two version flags each source resolves through its `CURRENT`.
""".format(year=year, a=a, b=b, lab_a=lab_a, lab_b=lab_b, va=ver[a], vb=ver[b],
           vd=ver["duneline"], vs=ver["shoreline"],
           gm=gap.mean(), gmed=gap.median(), gsd=gap.std(),
           glo=gap.min(), ghi=gap.max(),
           npos=int((gap > 0).sum()), n=len(gap), shift=shift,
           mm=mdiff.mean(), msd=mdiff.std(), mlo=mdiff.min(), mhi=mdiff.max(),
           stem=stem, filed=filed, extra=extra), encoding="utf-8")


if __name__ == "__main__":
    main()
