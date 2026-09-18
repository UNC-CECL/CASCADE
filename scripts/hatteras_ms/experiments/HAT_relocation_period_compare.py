#!/usr/bin/env python3
"""One relocation event, two hindcast windows: does the start year change
whether CASCADE reproduces it?

THE QUESTION
    The 1999 NC-12 relocation (GIS 9-14) sits inside BOTH the 1984-2004 and
    the 1996-2010 hindcast windows. Each window has its own emergent-vs-
    prescribed comparison (HAT_relocation_comparison.py, one set per preset),
    scored on its own terms. This script reads those two sets side by side
    for the domains ONE event moved and asks what the start year changed:

      * how much dune retreat each window accumulates at the road before
        the event year -- 15 model years from a 1984 start, 3 from 1996;
      * whether the free-running arm fires at all, and when, in each;
      * how each window's modelled 2004 road position compares with the one
        surveyed position, which is the END of one window and the MIDDLE of
        the other.

    It re-scores nothing. Every number is read from the per-period tables,
    so a disagreement between this report and a per-period report is a
    stale set, not a second opinion.

WHAT THE TWO WINDOWS SHARE, AND WHAT THEY DO NOT
    Same topography product (1984-start, one dune-topo version, read from the
    sets and required to agree), same road line (1978 for both starts), same
    code, same relocation target. They differ in the start year, the storm
    series, the offset survey (1984 line vs the 1997 line), and the setback
    file: a 1996 start reads the 1984 setbacks with the 1989 event already
    applied, which does not touch GIS 9-14. So at the event domains the two
    windows start the road in the SAME place and differ only in how many
    years of modelled retreat precede 1999.

USAGE
    python scripts/hatteras_ms/experiments/HAT_relocation_period_compare.py
    python scripts/hatteras_ms/experiments/HAT_relocation_period_compare.py --presets zeroBE edgeBE --version v2

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import argparse
import datetime
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
PROJECT_BASE_DIR = next(_p for _p in _HERE.parents if (_p / "pyproject.toml").exists())
SCRIPTS_DIR = PROJECT_BASE_DIR / "scripts"
for _path in (SCRIPTS_DIR, _HERE.parent):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from hatteras_site_config import HATTERAS_PERIODS, HATTERAS_ROAD_EVENTS   # noqa: E402
from cascade_pipeline.roadway import RelocationEvent                      # noqa: E402
import hat_figure_style as style                                          # noqa: E402

COMPARISONS = PROJECT_BASE_DIR / "output" / "comparisons"
CHECK_YEAR = 2004                 # the one surveyed road position
TOLERANCE_YEARS = (2, 5)
DEFAULT_PERIODS = (1984, 1996)
DEFAULT_PRESETS = ("zeroBE", "edgeBE")
DEFAULT_EVENT = 1999


# =============================================================================
# Reading the per-period sets
# =============================================================================

def set_dir(start_year, version, preset):
    """The per-period comparison folder this reads: relocation_<s>_<e>/<version>/<preset>/."""
    end = HATTERAS_PERIODS[start_year]["end_year"]
    return COMPARISONS / f"relocation_{start_year}_{end}" / version / preset


def event_domains(event_year):
    """The GIS domains one relocation event moved, from HATTERAS_ROAD_EVENTS."""
    for ev in HATTERAS_ROAD_EVENTS:
        if isinstance(ev, RelocationEvent) and ev.enabled and ev.year == event_year:
            return sorted(ev.displacement_m)
    raise SystemExit(f"no enabled relocation event in {event_year}")


def read_set(folder):
    """The tables and the provenance header of one per-period set."""
    tables = folder / "tables"
    need = ["first_relocation_year.csv", "setback_by_year.csv",
            "setback_summary.csv", "road_outcomes.csv", "confusion.csv"]
    missing = [n for n in need if not (tables / n).exists()]
    if missing:
        raise SystemExit(f"{folder}: missing {missing}; run "
                         "HAT_relocation_comparison.py for this period first")
    report = folder / "report.txt"
    header = []
    if report.exists():
        for line in report.read_text(encoding="utf-8").splitlines()[1:12]:
            if line.startswith("=") or not line.strip():
                break
            header.append(line)
    return {
        "folder": folder,
        "header": header,
        "first": pd.read_csv(tables / "first_relocation_year.csv"),
        "by_year": pd.read_csv(tables / "setback_by_year.csv"),
        "summary": pd.read_csv(tables / "setback_summary.csv"),
        "outcomes": pd.read_csv(tables / "road_outcomes.csv"),
        "confusion": pd.read_csv(tables / "confusion.csv"),
    }


def _topo_version_of(header):
    """'v2' from a report header's arm lines, or 'v?'."""
    for line in header:
        if "topo " in line:
            tok = line.split("topo ", 1)[1].split()[0]      # 1984-start/v2
            return tok.split("/")[-1]
    return "v?"


# =============================================================================
# The per-domain comparison
# =============================================================================

def domain_table(sets, domains, event_year):
    """One row per (period, domain): retreat before the event, the free arm's
    answer, and the position check. Everything read, nothing re-scored."""
    rows = []
    for start, s in sets.items():
        by = s["by_year"]
        first = s["first"].set_index("gis")
        summ = s["summary"].set_index("gis")
        for gis in domains:
            d = by[by["gis"] == gis].sort_values("year")
            free = d.set_index("year")["setback_free_m"]
            pre = event_year - 1
            start_m = float(free.iloc[0]) if len(free) else np.nan
            pre_m = float(free.get(pre, np.nan))
            f = first.loc[gis] if gis in first.index else None
            m = summ.loc[gis] if gis in summ.index else None
            rows.append(dict(
                period=f"{start}-{HATTERAS_PERIODS[start]['end_year']}",
                start_year=start,
                gis=gis,
                years_before_event=event_year - start,
                start_setback_m=start_m,
                setback_before_event_m=pre_m,
                retreat_before_event_m=start_m - pre_m,
                modelled_first_year=(None if f is None or pd.isna(f["modelled_first_year"])
                                     else int(f["modelled_first_year"])),
                error_years=(None if f is None or pd.isna(f["error_years"])
                             else int(f["error_years"])),
                outcome=None if f is None else f["outcome"],
                min_setback_m=None if f is None else f["min_setback_m"],
                migration_needed_m=None if f is None else f["migration_needed_m"],
                free_at_2004_m=None if m is None else m.get("free_at_check_m"),
                prescribed_at_2004_m=None if m is None else m.get("prescribed_at_check_m"),
                measured_2004_m=None if m is None else m.get("measured_2004_m"),
            ))
    df = pd.DataFrame(rows)
    for arm in ("free", "prescribed"):
        df[f"{arm}_2004_error_m"] = df[f"{arm}_at_2004_m"] - df["measured_2004_m"]
    return df


def event_recall(dom, tolerance):
    """Hits among the event domains, from the FIRST modelled relocation year.
    The per-period confusion.csv counts every event in the window and any
    relocation year; this restricts to one event and uses the first firing,
    which is the model's answer to 'when did the dune reach the road'."""
    out = []
    for period, g in dom.groupby("period", sort=False):
        err = g["error_years"].astype(float)
        hits = int((err.abs() <= tolerance).sum())
        out.append(dict(period=period, tolerance_years=tolerance,
                        event_domains=len(g), hits=hits,
                        recall=hits / len(g) if len(g) else np.nan,
                        hit_domains=" ".join(str(x) for x in g.loc[err.abs() <= tolerance, "gis"])))
    return pd.DataFrame(out)


def outcomes_at(sets, domains):
    rows = []
    for start, s in sets.items():
        o = s["outcomes"]
        o = o[o["gis"].isin(domains)]
        for arm, g in o.groupby("arm"):
            rows.append(dict(
                period=f"{start}-{HATTERAS_PERIODS[start]['end_year']}", arm=arm,
                domains=len(g), drowned=int(g["drowned"].sum()),
                relocation_blocked=int(g["relocation_blocked"].sum()),
                relocations=int(g["relocations"].sum()),
                overwash_removed_m3=float(g["overwash_removed_m3"].sum()),
                dunes_rebuilt=int(g["dunes_rebuilt"].sum())))
    return pd.DataFrame(rows)


# =============================================================================
# Figure: the setback trajectories, both windows on one axis per domain
# =============================================================================

def trajectory_figure(sets, domains, event_year, preset, out_path):
    """One panel per event domain. Each window is a colour (the vintage pair:
    the earlier start in the RdBu red, the later in the blue); the free arm
    is solid, the prescribed arm dashed; the event year is a vertical rule;
    the surveyed 2004 position is a black marker. A line stops where that
    arm stopped managing the road."""
    style.apply_style()
    C = style.C
    colours = dict(zip(sorted(sets), (C["EARLY"], C["LATE"])))
    ncol = 2
    nrow = int(np.ceil(len(domains) / ncol))
    fig, axes = plt_subplots(nrow, ncol, style.figsize("double", aspect=0.42 * nrow))
    years_all = []
    for i, gis in enumerate(domains):
        ax = axes[i]
        measured = None
        for start, s in sets.items():
            by = s["by_year"]
            d = by[by["gis"] == gis].sort_values("year")
            years_all += d["year"].tolist()
            col = colours[start]
            label_p = f"{start}-{HATTERAS_PERIODS[start]['end_year']}"
            fm = d[d["managed_free"]]
            pm = d[d["managed_prescribed"]]
            ax.plot(fm["year"], fm["setback_free_m"], color=col, lw=1.4,
                    label=f"{label_p}, emergent")
            ax.plot(pm["year"], pm["setback_prescribed_m"], color=col, lw=1.1,
                    ls="--", label=f"{label_p}, prescribed")
            summ = s["summary"].set_index("gis")
            if gis in summ.index and pd.notna(summ.loc[gis].get("measured_2004_m")):
                measured = float(summ.loc[gis]["measured_2004_m"])
        ax.axvline(event_year, color=C["INK_MUTED"], lw=0.7, ls=":", zorder=0)
        ax.axhline(0, color=C["INK"], lw=0.5, zorder=0)
        if measured is not None:
            ax.plot([CHECK_YEAR], [measured], marker="x", ms=6, mew=1.2,
                    color=C["ROAD"], ls="none", label="surveyed 2004 position",
                    zorder=5)
        style._title(ax, i, f"GIS {gis}")
        ax.grid(True, axis="y")
        ax.set_axisbelow(True)
    for ax in axes[len(domains):]:
        ax.set_visible(False)
    lo, hi = min(years_all), max(years_all)
    for ax in axes[:len(domains)]:
        ax.set_xlim(lo - 0.5, hi + 0.5)
    for r in range(nrow):
        axes[r * ncol].set_ylabel("NC-12 setback (m)")
    for ax in axes[max(0, len(domains) - ncol):len(domains)]:
        ax.set_xlabel("Year")
    handles, labels = axes[0].get_legend_handles_labels()
    seen, h2, l2 = set(), [], []
    for h, l in zip(handles, labels):
        if l not in seen:
            seen.add(l); h2.append(h); l2.append(l)
    fig.legend(h2, l2, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    return style.save(fig, out_path, close=True)


def plt_subplots(nrow, ncol, size):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(nrow, ncol, figsize=size, sharey=True)
    return fig, np.atleast_1d(axes).ravel()


# =============================================================================
# Report
# =============================================================================

def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--periods", type=int, nargs=2, default=DEFAULT_PERIODS,
                    help="two hindcast start years (default 1984 1996)")
    ap.add_argument("--presets", nargs="+", default=list(DEFAULT_PRESETS))
    ap.add_argument("--version", default="v2",
                    help="dune-topo version folder both sets were made on")
    ap.add_argument("--event", type=int, default=DEFAULT_EVENT)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    domains = event_domains(args.event)
    out_dir = (Path(args.out) if args.out else
               COMPARISONS / "relocation_periods" / f"{args.event}_event" / args.version)
    tables_dir = out_dir / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    lines = []

    def say(text=""):
        print(text)
        lines.append(text)

    say("=" * 74)
    say(f"generated   {datetime.datetime.now():%Y-%m-%d %H:%M:%S} by {_HERE.name}")
    say(f"event       {args.event}  GIS {domains}")
    say(f"periods     {' and '.join(f'{p}-{HATTERAS_PERIODS[p]['end_year']}' for p in args.periods)}")
    say(f"version     {args.version}")
    say("Read from the per-period sets named below; nothing re-scored.")
    say("=" * 74)

    captions = []
    for preset in args.presets:
        sets = {p: read_set(set_dir(p, args.version, preset)) for p in args.periods}
        versions = {p: _topo_version_of(s["header"]) for p, s in sets.items()}
        say()
        say("#" * 74)
        say(f"PRESET {preset}")
        say("#" * 74)
        for p, s in sets.items():
            say(f"  {p}: {s['folder'].relative_to(PROJECT_BASE_DIR)}")
            for h in s["header"]:
                say(f"      {h}")
        if len(set(versions.values())) > 1:
            say(f"  !! the two sets are on different dune-topo versions: {versions}")

        dom = domain_table(sets, domains, args.event)
        dom.to_csv(tables_dir / f"event_domains_{preset}.csv", index=False)

        say()
        say("-" * 74)
        say("1. RETREAT ACCUMULATED BEFORE THE EVENT, free-running arm")
        say("-" * 74)
        say("  start_setback: the road at year 0 (same file for both starts at these")
        say(f"  domains). setback_before_event: the free arm in {args.event - 1}.")
        cols = ["period", "gis", "years_before_event", "start_setback_m",
                "setback_before_event_m", "retreat_before_event_m"]
        say(dom[cols].to_string(index=False))
        for period, g in dom.groupby("period", sort=False):
            say(f"  {period}: mean retreat before {args.event} "
                f"{g['retreat_before_event_m'].mean():+.0f} m over "
                f"{int(g['years_before_event'].iloc[0])} yr")

        say()
        say("-" * 74)
        say("2. DID THE FREE ARM FIRE, AND WHEN")
        say("-" * 74)
        cols = ["period", "gis", "modelled_first_year", "error_years", "outcome",
                "min_setback_m", "migration_needed_m"]
        say(dom[cols].to_string(index=False))
        rec = pd.concat([event_recall(dom, t) for t in TOLERANCE_YEARS])
        rec.to_csv(tables_dir / f"event_recall_{preset}.csv", index=False)
        say()
        say(f"  recall on the {len(domains)} event domains, from the first modelled year:")
        say(rec.to_string(index=False))
        for p, s in sets.items():
            c = s["confusion"]
            fp = c["false_positives"].iloc[0]
            ctrl = c["control_domains"].iloc[0]
            say(f"  {p}: false positives period-wide {fp}/{ctrl} control domains")

        say()
        say("-" * 74)
        say(f"3. POSITION CHECK AT {CHECK_YEAR}")
        say("-" * 74)
        for p in args.periods:
            end = HATTERAS_PERIODS[p]["end_year"]
            say(f"  {p}-{end}: {CHECK_YEAR} is year {CHECK_YEAR - p} of {end - p}"
                + ("  (the end)" if CHECK_YEAR == end else "  (mid-window)"))
        cols = ["period", "gis", "free_at_2004_m", "prescribed_at_2004_m",
                "measured_2004_m", "free_2004_error_m", "prescribed_2004_error_m"]
        say(dom[cols].to_string(index=False))
        for period, g in dom.groupby("period", sort=False):
            say(f"  {period}: mean |error| free {g['free_2004_error_m'].abs().mean():.0f} m, "
                f"prescribed {g['prescribed_2004_error_m'].abs().mean():.0f} m")

        say()
        say("-" * 74)
        say("4. ROAD OUTCOMES AT THE EVENT DOMAINS")
        say("-" * 74)
        oc = outcomes_at(sets, domains)
        oc.to_csv(tables_dir / f"event_outcomes_{preset}.csv", index=False)
        say(oc.to_string(index=False))

        fig_path = out_dir / f"setback_trajectories_{args.event}_{preset}"
        written = trajectory_figure(sets, domains, args.event, preset, fig_path)
        say(f"\n  figure -> {written[0].name}")
        captions.append(
            f"**{fig_path.name}.png** -- NC-12 setback behind the dune line at each "
            f"domain the {args.event} relocation moved (GIS {domains[0]}-{domains[-1]}), "
            f"under the {preset} source/sink preset, for the two hindcast windows that "
            f"contain the event. Red is the {args.periods[0]} start, blue the "
            f"{args.periods[1]} start; solid is the emergent arm (the roadway module "
            f"deciding alone), dashed the prescribed arm (the measured {args.event} "
            f"displacement applied). The dotted rule is the event year, the cross the "
            f"surveyed {CHECK_YEAR} road position (RoadOffset_2004, the same "
            f"observation for both windows). A line ends where that arm stopped "
            f"managing the road. Setbacks move in whole 10 m cells; the emergent "
            f"trigger fires when the setback goes below zero.")

    report = out_dir / "report.txt"
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    (out_dir / "CAPTIONS.md").write_text(
        f"# Captions -- {args.event} relocation across hindcast windows\n\n"
        + "\n\n".join(captions) + "\n", encoding="utf-8")
    print(f"\nreport -> {report}")


if __name__ == "__main__":
    main()
