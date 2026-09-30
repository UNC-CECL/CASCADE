#!/usr/bin/env python3
"""
Redraw a finished run's figures without re-running the model, from its saved shoreline matrix.

    python rerender_run_figures.py --dry-run
    python rerender_run_figures.py --arm matrix/1996_2010/edgeBE [--gifs] [--open-npz]
    python rerender_run_figures.py --run-dir output/raw_runs/.../HAT_... [--ylim=-10,10]

Overwrites only image files it can rebuild exactly; never touches the .npz,
the matrix, the tables or run_index.csv. GIFs with relocation markers need
--open-npz. Details: scripts/figure_making/model_output/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""
from __future__ import annotations

import argparse
import dataclasses
import csv
import fnmatch
import json
import os
import sys
import time
import traceback
from pathlib import Path

import numpy as np


# The project root: the first parent holding data/hatteras_init
def _find_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit("could not locate the project root")


REPO = _find_root(Path(__file__).resolve())
sys.path.insert(0, str(REPO / "scripts"))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from cascade_pipeline.coastsat_lowess import (  # noqa: E402
    CoastSatDataset, LowessConfig, build_coastsat_series, scale_coastsat_series)
from cascade_pipeline.hindcast import build_shoreline_target  # noqa: E402
from cascade_pipeline.run_info import RunInfo  # noqa: E402
from cascade_pipeline.run_layout import ANIMATIONS, resolve  # noqa: E402
from cascade_pipeline.shoreline import compute_change_rate, compute_lrr  # noqa: E402
from cascade_pipeline.plotting.rate_comparison import (  # noqa: E402
    DEFAULT_RATE_COMPARISON, plot_annotated_rate_comparison,
    plot_rate_comparison)
from cascade_pipeline.plotting.shoreline_gif import (  # noqa: E402
    GifConfig, make_all_shoreline_gifs)
from site_layer.hatteras_site_config import (  # noqa: E402
    HATTERAS_ANNOTATIONS, HATTERAS_DOMAINS)

# --- CONFIG ------------------------------------------------------------------
RAW_RUNS = REPO / "output" / "raw_runs"
# Resolved through hat_observed_rates.py (2026-09-18), not typed.
from site_layer.hat_observed_rates import COASTSAT_LRR_ROOT as COASTSAT_BASE_DIR  # noqa: E402
from site_layer.hat_topo_version import RAW_OFFSET_DIR  # noqa: E402

# Must match section 8/9 of the hindcast runner; check_conventions asserts it
LOWESS_CONFIG = LowessConfig(window_domains=(7,), skip_southern_domains=10)
RATE_ESTIMATOR = "lrr"
FLIP_SIGN_MODEL = True
PLOT_REAL_DOMAINS_ONLY = True

COASTSAT_DATASETS = [
    CoastSatDataset(
        label="CoastSat LRR (1984-2004)", period_start=1984,
        csv_path=str(COASTSAT_BASE_DIR / "1984_2004" / "transect_lrr_full.csv")),
    CoastSatDataset(
        label="CoastSat LRR (2004-2024)", period_start=2004,
        csv_path=str(COASTSAT_BASE_DIR / "2004_2024" / "transect_lrr_full.csv")),
    # The 1996 and 2010 windows, as the runner lists them
    CoastSatDataset(
        label="CoastSat LRR (1996-2010)", period_start=1996,
        csv_path=str(COASTSAT_BASE_DIR / "1996_2010" / "transect_lrr_full.csv")),
    CoastSatDataset(
        label="CoastSat LRR (2010-2024)", period_start=2010,
        csv_path=str(COASTSAT_BASE_DIR / "2010_2024" / "transect_lrr_full.csv")),
]

# The full-record rate the projected position change is built from
LONG_TERM_DATASET = CoastSatDataset(
    label="CoastSat LRR (1996-2024)", period_start=1996,
    csv_path=str(COASTSAT_BASE_DIR / "1996_2024" / "transect_lrr_full.csv"))
LONG_TERM_WINDOW = (1996, 2024)
POSITION_DIR = "position_change"


# The two observed references a run's position change is drawn against
def position_change_jobs(run, cs_cache):
    span = run.end_year - run.start_year
    win = f"{run.start_year}–{run.end_year}"
    lt = f"{LONG_TERM_WINDOW[0]}–{LONG_TERM_WINDOW[1]}"
    widest = max(LOWESS_CONFIG.window_domains)
    if "long_term" not in cs_cache:
        cs_cache["long_term"] = build_coastsat_series(
            [LONG_TERM_DATASET], active_period_start=LONG_TERM_WINDOW[0],
            lowess_config=LOWESS_CONFIG, domains=HATTERAS_DOMAINS)
    return [
        ("total",
         scale_coastsat_series(cs_cache[run.start_year], span),
         f"CoastSat total change (LRR {win} × {span} yr, {widest}-domain LOWESS)",
         f"total change, the per-transect LRR fitted on {win} multiplied by "
         f"{span} yr"),
        ("projected",
         scale_coastsat_series(cs_cache["long_term"], span, active=True),
         f"CoastSat projected change (LRR {lt} × {span} yr, {widest}-domain LOWESS)",
         f"projected change, the per-transect LRR fitted on the full "
         f"{lt} record multiplied by {span} yr"),
    ]


GIF_JOBS = [
    dict(range="real", mode="displacement"),
    dict(range="real", mode="position"),
    dict(range="groin", mode="position", pad=9),
    dict(range="groin", mode="difference", pad=9),
]
# -----------------------------------------------------------------------------


# Stop if the runner's figure conventions no longer match the constants here
def check_conventions() -> None:
    src = (REPO / "scripts" / "hatteras_ms" / "HAT_hindcast_1984_2024.py").read_text(
        encoding="utf-8", errors="replace")
    want = {
        'RATE_ESTIMATOR = "lrr"': RATE_ESTIMATOR == "lrr",
        "window_domains=(7,)": LOWESS_CONFIG.window_domains == (7,),
        "skip_southern_domains=10": LOWESS_CONFIG.skip_southern_domains == 10,
        "PLOT_REAL_DOMAINS_ONLY = True": PLOT_REAL_DOMAINS_ONLY is True,
        "FLIP_SIGN_MODEL = True": FLIP_SIGN_MODEL is True,
    }
    missing = [k for k, ok in want.items() if not (ok and k in src)]
    if missing:
        raise SystemExit(
            "the hindcast's figure conventions no longer match this script:\n  "
            + "\n  ".join(missing)
            + "\nUpdate the constants at the top of this file to match "
              "HAT_hindcast_1984_2024.py before re-rendering.")


# Run folders to act on: those holding a *_run_metadata.json
def find_runs(args) -> list[Path]:
    if args.run_dir:
        roots = [Path(args.run_dir).resolve()]
    else:
        roots = sorted({p.parent for p in RAW_RUNS.rglob("*_run_metadata.json")})
    out = []
    for d in roots:
        rel = d.relative_to(RAW_RUNS).as_posix() if RAW_RUNS in d.parents else d.as_posix()
        if args.arm and not rel.startswith(args.arm.strip("/")):
            continue
        if args.match and not fnmatch.fnmatch(d.name, args.match):
            continue
        out.append(d)
    return out


# (run_name, metadata, shoreline matrix in m) for one run, or None if unusable
def load_run(run_dir: Path):
    meta_hits = sorted(run_dir.glob("*_run_metadata.json"))
    if not meta_hits:
        return None, None, None, "no run metadata"
    run_name = meta_hits[0].name[: -len("_run_metadata.json")]
    with open(meta_hits[0], encoding="utf-8") as fh:
        meta = json.load(fh)
    npy = resolve(run_dir, "matrix", run_name)
    if not npy.is_file():
        return run_name, meta, None, "no shoreline matrix"
    shoreline_m = np.load(npy)
    if shoreline_m.ndim != 2 or shoreline_m.shape[0] < 2:
        return run_name, meta, None, f"matrix shape {shoreline_m.shape}"
    return run_name, meta, shoreline_m, None


# The RunInfo the plotting package expects, from a run's metadata
def run_info_from(meta: dict, run_name: str, run_dir: Path) -> RunInfo:
    period = meta.get("period", {})
    wave = meta.get("wave climate", {})
    src = meta.get("source/sink", {})
    return RunInfo(
        run_name=run_name,
        run_dir=str(run_dir),
        start_year=int(period["start_year"]),
        end_year=int(period["end_year"]),
        Hs=float(wave["wave_height_m"]) if wave.get("wave_height_m") is not None else None,
        flip_sign_model=FLIP_SIGN_MODEL,
        background_erosion_on=bool(src.get("background_erosion_on", True)),
        wave_climate=meta.get("scenario", {}).get("wave climate"),
    )


# Do this run's GIFs carry road-relocation markers? (from its road table)
def has_relocation_markers(run_dir: Path, run_name: str) -> bool:
    path = resolve(run_dir, "road_csv", run_name)
    if not path.is_file():
        return False
    try:
        with open(path, newline="", encoding="utf-8") as fh:
            return any(int(float(r.get("relocations") or 0)) > 0
                       for r in csv.DictReader(fh))
    except (ValueError, KeyError):
        return False


# Per-year relocation flags, by opening the archive (slow)
def relocations_from_npz(run_dir: Path, run_name: str, n_states: int):
    npz = resolve(run_dir, "archive", run_name)
    if not npz.is_file():
        return None
    cascade = np.load(npz, allow_pickle=True)["cascade"][0]
    roadways = getattr(cascade, "_roadways", None)
    if not roadways:
        return None
    events = np.zeros((n_states, HATTERAS_DOMAINS.total_domains), dtype=bool)
    for pad, roadway in enumerate(roadways):
        raw = getattr(roadway, "_road_relocated_TS", None)
        series = (np.asarray([], dtype=float) if raw is None
                  else np.asarray(raw, dtype=float))
        if series.size:
            n = min(series.size, events.shape[0])
            events[:n, pad] = series[:n] > 0
    return events


# Redraw one run's rate figures, position-change figures and (optionally) GIFs
def rerender(run_dir: Path, args, cs_cache: dict) -> dict:
    run_name, meta, shoreline_m, why = load_run(run_dir)
    result = {"dir": run_dir, "run": run_name, "figures": 0, "gifs": 0,
              "skipped": None, "error": None}
    if shoreline_m is None:
        result["skipped"] = why
        return result

    run = run_info_from(meta, run_name, run_dir)
    run_years = int(meta["period"].get("run_years", shoreline_m.shape[0] - 1))

    # The same two calls the run made, on the same saved matrix.
    if RATE_ESTIMATOR == "lrr":
        plotted_rate, _ = compute_lrr(shoreline_m, span_years=run_years,
                                      flip_sign=FLIP_SIGN_MODEL)
    else:
        plotted_rate = compute_change_rate(shoreline_m, span_years=run_years,
                                           flip_sign=FLIP_SIGN_MODEL)

    # CoastSat is the same for every run of a given start year; build once.
    key = run.start_year
    if key not in cs_cache:
        cs_cache[key] = build_coastsat_series(
            COASTSAT_DATASETS, active_period_start=key,
            lowess_config=LOWESS_CONFIG, domains=HATTERAS_DOMAINS)
    cs_series = cs_cache[key]

    config = dataclasses.replace(DEFAULT_RATE_COMPARISON, ylim=args.ylim,
                                 ylim_real=args.ylim_real)
    if args.lowess_only:
        config = dataclasses.replace(config, plot_domain_means=False,
                                     raw_lrr_southern_only=True,
                                     overlay_raw_southern_only=True,
                                     south_mean_line=True, show_features=True,
                                     model_label_plain=True,
                                     show_wave_climate=True,
                                     publication_text=True)
    fig_kwargs = dict(domains=HATTERAS_DOMAINS, annotations=HATTERAS_ANNOTATIONS,
                      lowess_config=LOWESS_CONFIG, config=config)
    # Resolved, not joined: overwrite the figure wherever the run keeps it
    rate_fig_kind = "figure_rate" if PLOT_REAL_DOMAINS_ONLY else "figure_rate_buffers"
    rate_png = resolve(run_dir, rate_fig_kind, run_name)
    buffers_png = resolve(run_dir, "figure_rate_buffers", run_name)

    if args.dry_run:
        result["figures"] = 2
        result["gifs"] = len(GIF_JOBS) if args.gifs else 0
        return result

    for _png in (rate_png, buffers_png):
        _png.parent.mkdir(parents=True, exist_ok=True)

    plot_rate_comparison(
        plotted_rate, cs_series, run,
        real_domains_only=PLOT_REAL_DOMAINS_ONLY, estimator=RATE_ESTIMATOR,
        save_path=str(rate_png),
        show=False, **fig_kwargs)
    plot_annotated_rate_comparison(
        plotted_rate, cs_series, run, estimator=RATE_ESTIMATOR,
        save_path=str(buffers_png),
        show=False, **fig_kwargs)
    result["figures"] = 2

    if args.position_change:
        # End minus start, in metres: compute_change_rate over a span of 1.
        position_m = compute_change_rate(shoreline_m, span_years=1,
                                         flip_sign=FLIP_SIGN_MODEL)
        for ref, series, obs_label, obs_text in position_change_jobs(run, cs_cache):
            pconfig = dataclasses.replace(
                config, quantity="position", observed_label=obs_label,
                observed_description=obs_text, ylim=None, ylim_real=None)
            pdir = rate_png.parent / POSITION_DIR / ref
            pdir.mkdir(parents=True, exist_ok=True)
            pkw = dict(fig_kwargs, config=pconfig)
            plot_rate_comparison(
                position_m, series, run,
                real_domains_only=PLOT_REAL_DOMAINS_ONLY, estimator=None,
                save_path=str(pdir / "shoreline_position_change.png"),
                show=False, **pkw)
            plot_annotated_rate_comparison(
                position_m, series, run, estimator=None,
                save_path=str(pdir / "shoreline_position_change_with_buffers.png"),
                show=False, **pkw)
            plt.close("all")
            result["figures"] += 2

    if args.gifs:
        # Refresh only: a run without GIFs does not get new ones
        if not (any(run_dir.glob("*.gif"))
                or any((run_dir / ANIMATIONS).glob("*.gif"))):
            result["skipped"] = "gifs: this run has none to refresh"
            plt.close("all")
            return result
        markers = has_relocation_markers(run_dir, run_name)
        relocations = None
        if markers:
            if not args.open_npz:
                result["skipped"] = ("gifs: carries relocation markers; "
                                     "re-run with --open-npz to keep them")
                plt.close("all")
                return result
            relocations = relocations_from_npz(run_dir, run_name,
                                               shoreline_m.shape[0])
            if relocations is None:
                result["skipped"] = "gifs: markers expected but no archive to read"
                plt.close("all")
                return result

        target_m, _ = build_shoreline_target(
            shoreline_m[0], run.start_year, run.end_year, HATTERAS_DOMAINS,
            RAW_OFFSET_DIR)
        baseline = run_dir.parent / f"{run_name.replace('_groin', '_nogroin')}"
        baseline_npy = resolve(baseline, "matrix", baseline.name)
        gif_config = GifConfig(
            fps=3, year_stride=1, annotate=True, auto_open=False,
            keep_frames=False,
            save_matrix=False,          # the matrix on disk is the run's own
            ocean_at_bottom=True,
            baseline_label="no-groin baseline",
            target_label=f"observed {run.end_year} dune line")
        paths = make_all_shoreline_gifs(
            shoreline_m, run, GIF_JOBS,
            baseline_npy=str(baseline_npy) if baseline_npy.is_file() else None,
            target_m=target_m, relocations=relocations,
            domains=HATTERAS_DOMAINS, annotations=HATTERAS_ANNOTATIONS,
            gif_config=gif_config)
        result["gifs"] = len(paths or [])

    plt.close("all")
    return result


# Run: find the runs, check conventions, redraw each, report
def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__.splitlines()[1],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--arm", default=None,
                    help="only runs whose path under output/raw_runs starts "
                         "with this, e.g. matrix/1984_2004/calibBE or "
                         "sensitivity/waveHs")
    ap.add_argument("--match", default=None,
                    help="glob on the run folder name, e.g. '*calibBE*groin'")
    ap.add_argument("--run-dir", default=None, help="one run directory")
    ap.add_argument("--gifs", action="store_true",
                    help="also redraw the animations (slow)")
    ap.add_argument("--open-npz", action="store_true",
                    help="with --gifs, open the archive for runs whose GIFs "
                         "carry relocation markers (240 MB per run)")
    ap.add_argument("--dry-run", action="store_true",
                    help="list what would be redrawn; write nothing")
    ap.add_argument("--limit", type=int, default=0, help="stop after N runs")
    ap.add_argument("--ylim", default=None,
                    type=lambda s: tuple(float(v) for v in s.split(",")),
                    help="LOW,HIGH in m/yr: one y axis for every run redrawn, "
                         "so a set can be compared figure to figure (e.g. "
                         "--ylim=-10,10); default fits each run to its data")
    ap.add_argument("--ylim-real", default=None,
                    type=lambda s: tuple(float(v) for v in s.split(",")),
                    help="LOW,HIGH for the real-domains-only figure alone "
                         "(default: --ylim), e.g. --ylim-real=-7.5,7.5")
    ap.add_argument("--lowess-only", action="store_true",
                    help="the sensitivity-figure style: LOWESS curve, plus "
                         "transect dots and a dashed unsmoothed domain-mean "
                         "line over the southern domains only; shoals, piers "
                         "and groin drawn; no title -- one tag with the period "
                         "and wave settings, the rest in supporting/CAPTIONS.md")
    ap.add_argument("--position-change", action="store_true",
                    help="with --lowess-only: also draw the run's shoreline "
                         "position change (end minus start, m) against CoastSat "
                         "total change (the window's own LRR x span) and "
                         "projected change (the 1996-2024 LRR x span), into "
                         "figures/position_change/{total,projected}/")
    ap.add_argument("--traceback", action="store_true",
                    help="print a full traceback for a failed run")
    args = ap.parse_args()

    check_conventions()
    runs = find_runs(args)
    if args.limit:
        runs = runs[: args.limit]
    if not runs:
        raise SystemExit("no run directories matched")

    print(f"{'DRY RUN: ' if args.dry_run else ''}{len(runs)} run(s)"
          f"{' , animations included' if args.gifs else ''}\n")
    cs_cache: dict = {}
    started = time.time()
    done = skipped = failed = 0
    for i, d in enumerate(runs, 1):
        rel = d.relative_to(RAW_RUNS).as_posix() if RAW_RUNS in d.parents else d.name
        try:
            r = rerender(d, args, cs_cache)
        except Exception as exc:                      # noqa: BLE001
            failed += 1
            print(f"[{i:3}/{len(runs)}] FAILED  {rel}\n          "
                  f"{type(exc).__name__}: {exc}")
            if args.traceback:
                traceback.print_exc()
            continue
        if r["skipped"] and not r["figures"]:
            skipped += 1
            print(f"[{i:3}/{len(runs)}] skipped {rel}  ({r['skipped']})")
        else:
            done += 1
            note = f"  ({r['skipped']})" if r["skipped"] else ""
            print(f"[{i:3}/{len(runs)}] {r['figures']} fig"
                  f"{f', {r_gifs} gif' if (r_gifs := r['gifs']) else ''}"
                  f"  {rel}{note}")

    mins = (time.time() - started) / 60.0
    print(f"\n{done} redrawn, {skipped} skipped, {failed} failed "
          f"in {mins:.1f} min")
    if not args.dry_run:
        print("the .npz, .npy, CSV and TXT files and run_index.csv were not "
              "touched")


if __name__ == "__main__":
    main()
