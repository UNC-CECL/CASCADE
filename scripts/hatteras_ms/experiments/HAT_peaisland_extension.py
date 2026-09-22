r"""
HAT_peaisland_extension.py -- does the buffer's orientation matter?
==============================================================================
THE QUESTION (Hannah, 2026-09-16). The hindcast models GIS 1-90 and pads
each end with 15 invented domains that extrapolate the local shoreline
slope, then bridge back to close BRIE's periodic ring; edgeBE pins GIS 1
and 90 to their observed rates with boundary source/sink terms that take up
whatever the buffer gets wrong (+32.2 and +10.0 m/yr in 1996-2010). What if the
buffer carried the REAL coast instead -- Pea Island north to GIS 115 --
with the ends re-solved there? (A one-domain southern extension was run and
removed the same evening: no domain polygon lies south of GIS 1.)

THE DESIGN
    geometries   n115 (GIS 1-115) against base
    offset modes asrun (the compressed planform every calibrated run uses)
                 and detrended (the planform at full strength; no calibrated
                 baseline exists, so base-detrended is solved here too)
    period       1996-2010, full_management, no groin, relocations off, Hs 2.5
    topography   1984-start CURRENT for GIS 1-90; the buffer profile beyond
    management   none on the extension (no road, fills, relocation, BE)
    stage 0      zeroBE on every member: the orientation effect with no
                 boundary term anywhere, against the zeroBE matrix run
    stage 1      the end domains solved by Newton steps (edgeBE, one probe
                 per step through HAT_BE_OVERRIDE), against the edgeBE
                 matrix run
    score        interior RMSE on GIS 2-89 against the SURVEYED target (the
                 runner's rmse_interior_m_yr, identical in meaning for every
                 geometry), the solved end values, and the rates on GIS 80-90

WHERE THINGS ARE
    inputs   2-brie-offset/1996/ext/<geometry>/          the offsets
             5-scr/3-rates/coastsat/lrr/1996_2010/ext/            the targets
    runs     output/raw_runs/experiments/2026-09-16-peaisland-ext/<member>/
             one member per <geometry>-<mode>: its 1996_2010/zeroBE/ run is
             stage 0, its step<k>/ folders are the Newton probes, and SOLVED
             names the step that stands as the solved run
    logs     output/raw_runs/experiments/2026-09-16-peaisland-ext/logs/<member>/
    answer   RESULTS.md and figures/ beside NOTE.md in that folder

USAGE
    python HAT_peaisland_extension.py stage0                  # the 5 zeroBE runs
    python HAT_peaisland_extension.py check                   # base geometry, edgeBE:
                                                              # must reproduce the matrix row
    python HAT_peaisland_extension.py probe --member n115-asrun --step 1 \
        --override "1=32.2,115=12.0"                          # one Newton probe
    python HAT_peaisland_extension.py next --member n115-asrun  # the next probe, from
                                                              # the runs so far
    python HAT_peaisland_extension.py score                   # RESULTS.md + figures

Each run is the ordinary hindcast runner driven through the environment,
exactly as HAT_run_all.py drives the matrix; nothing here reimplements a
run. Runs are never overwritten (pass --overwrite to redo one).
==============================================================================
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

_HERE = Path(__file__).resolve()
PROJECT_ROOT = next(_p for _p in _HERE.parents
                    if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from site_layer.hat_extension_domains import GEOMETRIES, BASE_GEOMETRY  # noqa: E402

HINDCAST = PROJECT_ROOT / "scripts" / "hatteras_ms" / "HAT_hindcast_1984_2024.py"
SOLVE = (PROJECT_ROOT / "scripts" / "input_prep" / "7-source-sink" / "2-calibrate"
         / "be_edge_domain_solve.py")
RAW_RUNS = PROJECT_ROOT / "output" / "raw_runs"

TAG = "2026-09-16-peaisland-ext"
EXPERIMENT_DIR = RAW_RUNS / "experiments" / TAG
LOG_DIR = EXPERIMENT_DIR / "logs"
PERIOD = 1996
SCENARIO = "full_management"
RUN_TIMEOUT_S = 3600

# member -> (geometry, offset mode)
MEMBERS = {
    "n115-asrun": ("n115", "asrun"),
    "n115-detrended": ("n115", "detrended"),
    # the detrended planform has no calibrated 90-domain run to compare
    # against, so its baseline is solved here alongside the extensions
    "base-detrended": (BASE_GEOMETRY, "detrended"),
    # the regression check: base geometry, compressed, edgeBE -- must land
    # on the matrix row's numbers exactly
    "base-check": (BASE_GEOMETRY, "asrun"),
}


def run_env(member, preset, override="", overwrite=False):
    """The environment one run reads, built the way HAT_run_all builds it:
    every HAT_* variable named here, none inherited from the shell."""
    geometry, mode = MEMBERS[member]
    env = {k: v for k, v in os.environ.items() if not k.startswith("HAT_")}
    env.update({
        "HAT_IGNORE_SETTINGS": "1",
        "HAT_START_YEAR": str(PERIOD),
        "HAT_SOURCE_SINK_PRESET": preset,
        "HAT_SCENARIO": SCENARIO,
        "HAT_RELOCATIONS": "false",
        "HAT_GROIN_ENABLED": "false",
        "HAT_GEOMETRY": geometry,
        "HAT_OFFSET_MODE": mode,
        "HAT_RUN_KIND": "experiment",
        "HAT_OVERWRITE": "true" if overwrite else "false",
        "HAT_SAVE_MODEL_STATE": "false",
        "HAT_MAKE_GIFS": "false",
        "MPLBACKEND": "Agg",
        "PYTHONIOENCODING": "utf-8",
    })
    if override:
        env["HAT_BE_OVERRIDE"] = override
    return env


def launch(member, preset, step=None, override="", overwrite=False, dry_run=False):
    """One hindcast run, filed under experiments/<TAG>/<member>/ (stage 0)
    or experiments/<TAG>/<member>/step<k>/ (a Newton probe)."""
    tag_member = member if step is None else f"{member}/step{step}"
    env = run_env(member, preset, override, overwrite)
    env["HAT_RUN_TAG"] = f"{TAG}/{tag_member}"
    label = f"{tag_member} {preset}" + (f" override {override}" if override else "")
    if dry_run:
        print(f"would run: {label}")
        return True
    (LOG_DIR / member).mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / member / (f"step{step}_{preset}.log" if step else f"stage0_{preset}.log")
    print(f"running {label} ... ", end="", flush=True)
    t0 = time.perf_counter()
    proc = subprocess.run([sys.executable, str(HINDCAST)], env=env,
                          cwd=str(PROJECT_ROOT), capture_output=True, text=True,
                          encoding="utf-8", errors="replace", timeout=RUN_TIMEOUT_S)
    log_path.write_text((proc.stdout or "") + "\n--- STDERR ---\n" + (proc.stderr or ""),
                        encoding="utf-8")
    minutes = (time.perf_counter() - t0) / 60
    if proc.returncode != 0:
        tail = "\n".join((proc.stderr or "").strip().splitlines()[-6:])
        print(f"FAILED after {minutes:.1f} min (exit {proc.returncode}); log {log_path}\n{tail}")
        return False
    print(f"done in {minutes:.1f} min; log {log_path.name}")
    return True


def cmd_stage0(a):
    ok = True
    for member in (a.members or [m for m in MEMBERS if m != "base-check"]):
        ok &= launch(member, "zeroBE", overwrite=a.overwrite, dry_run=a.dry_run)
    return 0 if ok else 1


def cmd_check(a):
    return 0 if launch("base-check", "edgeBE", overwrite=a.overwrite, dry_run=a.dry_run) else 1


def cmd_probe(a):
    return 0 if launch(a.member, "edgeBE", step=a.step, override=a.override,
                       overwrite=a.overwrite, dry_run=a.dry_run) else 1


def solve_history(member):
    """(run_name, preset, tag) of the member's stage-0 run and every probe
    on disk, oldest first, read off the experiment folder."""
    from cascade_pipeline.run_registry import load_run_index
    index = load_run_index(RAW_RUNS / "run_index.csv")
    mine = index[(index["kind"] == "experiment")
                 & (index["tag"].str.startswith(f"{TAG}/{member}"))]
    rows = []
    for _, r in mine.iterrows():
        m, s = _member_step(r["tag"])
        if m != member:
            continue
        rows.append((s, r["run_name"], r["source_sink_preset"], r["tag"]))
    return [r[1:] for r in sorted(rows)]


def _member_step(tag):
    """(member, step) from a run's tag: <TAG>/<member> is stage 0,
    <TAG>/<member>/step<k> the k-th Newton probe. (None, None) otherwise."""
    parts = tag.split("/")
    if len(parts) == 2 and parts[1] in MEMBERS:
        return parts[1], 0
    if len(parts) == 3 and parts[1] in MEMBERS and parts[2].startswith("step"):
        return parts[1], int(parts[2][4:])
    return None, None


def mark_solved(member, step):
    """A SOLVED file in the member folder naming the step that is the
    solved run, for a reader who is not going to parse the index."""
    (EXPERIMENT_DIR / member).mkdir(parents=True, exist_ok=True)
    (EXPERIMENT_DIR / member / "SOLVED").write_text(f"step{step}\n", encoding="utf-8")


def next_probe(member, quiet=False):
    """The solve script over the member's runs so far, with the member's
    geometry in the environment. Returns (override or None, converged,
    text)."""
    import re
    history = solve_history(member)
    if not history:
        sys.exit(f"no runs on disk for {member}; run stage0 first")
    geometry, _ = MEMBERS[member]
    env = {k: v for k, v in os.environ.items() if not k.startswith("HAT_")}
    env["HAT_GEOMETRY"] = geometry
    env["PYTHONIOENCODING"] = "utf-8"
    # the zeroBE stage-0 run sits under zeroBE/, the probes under edgeBE/;
    # the solve script reads each run's preset off the index, so every run
    # is passed with its own tag and nothing else
    cmd = [sys.executable, str(SOLVE), "--period", str(PERIOD), "--kind", "experiment"]
    for run_name, _preset, tag in history:
        cmd += ["--run", run_name, "--tag", tag]
    proc = subprocess.run(cmd, env=env, cwd=str(PROJECT_ROOT), text=True,
                          encoding="utf-8", errors="replace", capture_output=True)
    text = (proc.stdout or "") + (proc.stderr or "")
    if not quiet:
        print("$ " + " ".join(cmd[2:]))
        print(text)
    if proc.returncode != 0:
        return None, False, text
    m = re.search(r'HAT_BE_OVERRIDE="([^"]+)"', text)
    # Per end: the solve script prints CONVERGED once the residual is inside
    # tolerance, and "no secant" once two probes imposed the same value --
    # which only happens after a converged end's step rounded to 0.0. Either
    # is done. Converged means every end is.
    blocks = re.split(r"^GIS ", text, flags=re.M)[1:]
    done = [("CONVERGED" in b) or ("no secant" in b) for b in blocks]
    converged = bool(done) and all(done)
    override = m.group(1) if m else None
    if override and not converged:
        # An end the script gave no step for (converged, or no secant)
        # keeps its last imposed value: the runner needs a rate at every end.
        given = {p.split("=")[0] for p in override.split(",")}
        last_imposed = {}
        for b in blocks:
            gis = b.split()[0]
            # a probe row is "<run name> <imposed> <model> <residual>"; the
            # "next probe" line also starts with HAT_ and has no numbers
            rows = [l.split() for l in b.splitlines()
                    if l.strip().startswith("HAT_") and len(l.split()) >= 4]
            if rows:
                last_imposed[gis] = rows[-1][-3]
        extra = [f"{g}={v}" for g, v in last_imposed.items()
                 if g not in given and float(v) != 0.0]
        if extra:
            override = ",".join([override] + extra)
    return override, converged, text


def cmd_next(a):
    override, converged, _ = next_probe(a.member)
    return 0 if override or converged else 1


# A probe whose |value| exceeds this is not a boundary term any more; the
# matrix values are tens of m/yr and Barrier3D's overwash router has been seen
# to die silently under runaway progradation. Stop and say so instead.
PROBE_CEILING_M_YR = 250.0


def cmd_solve(a):
    """Newton steps for one member until both ends converge or --max-steps
    is reached: next probe, run it, repeat."""
    history = solve_history(a.member)
    step = max((int(t.rsplit("step", 1)[-1]) for _, _, t in history if "-step" in t),
               default=0)
    for _ in range(a.max_steps):
        override, converged, text = next_probe(a.member, quiet=True)
        tail = "\n".join(l for l in text.splitlines()
                         if l.startswith("GIS") or "next" in l or "CONVERGED" in l
                         or "residual" in l and "run" not in l)
        print(f"[{a.member}] after step {step}:\n{tail}")
        if converged:
            print(f"[{a.member}] converged at step {step}")
            mark_solved(a.member, step)
            return 0
        if not override:
            print(f"[{a.member}] no next probe (see the solve output above)")
            print(text)
            return 1
        values = [abs(float(p.split("=")[1])) for p in override.split(",")]
        if max(values) > PROBE_CEILING_M_YR:
            print(f"[{a.member}] refusing probe {override}: beyond "
                  f"{PROBE_CEILING_M_YR:.0f} m/yr, not a boundary term any more")
            return 2
        step += 1
        if not launch(a.member, "edgeBE", step=step, override=override,
                      overwrite=a.overwrite, dry_run=a.dry_run):
            return 1
        if a.dry_run:
            return 0
    print(f"[{a.member}] stopped after {a.max_steps} step(s) without converging; "
          f"step {step} stands as the solved run until another is taken")
    mark_solved(a.member, step)
    return 3


# =============================================================================
# score: RESULTS.md and the figures
# =============================================================================

BASELINE_RUNS = {
    # (preset) -> the matrix run every compressed member is compared against
    "zeroBE": "HAT_1996_2010_zeroBE_road_bdm_nogroin",
    "edgeBE": "HAT_1996_2010_edgeBE_road_bdm_nogroin",
}
NEAR_END_GIS = (80, 90)     # the domains reported beside the interior score


def _rates(run_dir, run_name):
    import pandas as pd
    from cascade_pipeline.run_layout import resolve
    frame = pd.read_csv(resolve(run_dir, "rate_csv", run_name))
    return frame.set_index("gis_domain")["lrr_m_yr"]


def _target(extended=False):
    """The CoastSat target as the runner builds it: the surveyed GIS 1-90
    table (what the interior score uses), or the extension's table over
    GIS 0-115 (what an extended run's ends are solved against)."""
    from cascade_pipeline.coastsat_loess import (CoastSatDataset, LoessConfig,
                                                 build_coastsat_series)
    from cascade_pipeline.domains import DEFAULT_DOMAINS, DomainGeometry
    from cascade_pipeline.hindcast import build_target_table
    from site_layer.hat_observed_rates import lrr_csv, lrr_csv_ext
    loess = LoessConfig(window_domains=(10,), skip_southern_domains=10)
    if extended:
        first, last = GEOMETRIES["n115"]
        domains = DomainGeometry(num_real_domains=last - first + 1, first_gis_id=first)
        csv_path = lrr_csv_ext(PERIOD, 2010)
    else:
        domains, csv_path = DEFAULT_DOMAINS, lrr_csv(PERIOD, 2010)
    series = build_coastsat_series(
        [CoastSatDataset(label="extended" if extended else "surveyed",
                         period_start=PERIOD, csv_path=str(csv_path))],
        PERIOD, loess, domains=domains)
    table = build_target_table(series[0], loess, domains, 10)
    return table.set_index("gis_domain")["target_lrr_m_yr"]


def _surveyed_target():
    return _target(extended=False)


def collect():
    """Every run of the experiment plus the two matrix baselines, as rows:
    member, step, preset, index row, and the per-domain LRR."""
    from cascade_pipeline.run_registry import load_run_index, find_run_dir
    index = load_run_index(RAW_RUNS / "run_index.csv")
    rows = []
    mine = index[(index["kind"] == "experiment") & (index["tag"].str.startswith(TAG + "/"))]
    for _, r in mine.iterrows():
        member, step = _member_step(r["tag"])
        if member is None:
            continue
        run_dir = find_run_dir(RAW_RUNS, r["run_name"], (PERIOD, 2010),
                               r["source_sink_preset"], kind="experiment", tag=r["tag"])
        rows.append(dict(member=member, step=step, preset=r["source_sink_preset"],
                         geometry=MEMBERS[member][0], mode=MEMBERS[member][1],
                         row=r, lrr=_rates(run_dir, r["run_name"])))
    base = index[(index["kind"] == "matrix") & (index["status"] == "current")
                 & (index["start_year"].astype(str) == str(PERIOD))]
    for preset, name in BASELINE_RUNS.items():
        hit = base[base["run_name"] == name]
        if hit.empty:
            print(f"  no current matrix run {name}")
            continue
        r = hit.iloc[-1]
        run_dir = find_run_dir(RAW_RUNS, name, (PERIOD, 2010), preset)
        rows.append(dict(member="base-asrun", step=0, preset=preset, geometry="base",
                         mode="asrun", row=r, lrr=_rates(run_dir, name)))
    return rows


def final_runs(rows):
    """Per (member, preset), the last Newton step (edgeBE) or the stage-0 run."""
    out = {}
    for d in rows:
        key = (d["member"], d["preset"])
        if key not in out or d["step"] > out[key]["step"]:
            out[key] = d
    return out


def cmd_score(a):
    import numpy as np
    finals = final_runs(collect())
    target = _surveyed_target()
    lo, hi = NEAR_END_GIS
    near = list(range(lo, hi + 1))

    def num(row, key):
        v = row.get(key, "")
        return float(v) if v not in ("", None) else float("nan")

    lines = [f"# {TAG}: results", "",
             "Interior RMSE is GIS 2-89 against the surveyed CoastSat target in "
             "every geometry. End values are what the run imposed (m/yr). "
             f"\"near-end\" is the RMSE over GIS {lo}-{hi} against the same target. "
             "Stage 0 is zeroBE; the solved rows are the last Newton probe.", "",
             "| member | geometry | mode | preset | step | interior RMSE | interior bias "
             "| reach RMSE | near-end RMSE | south end | north end |",
             "|---|---|---|---|---|---|---|---|---|---|---|"]
    order = ["base-asrun", "base-check", "base-detrended", "n115-asrun", "n115-detrended"]
    for member in order:
        for preset in ("zeroBE", "edgeBE"):
            d = finals.get((member, preset))
            if d is None:
                continue
            r = d["row"]
            first, last = GEOMETRIES[d["geometry"]]
            resid = np.array([d["lrr"].get(g, np.nan) - target.get(g, np.nan) for g in near])
            near_rmse = float(np.sqrt(np.nanmean(resid ** 2)))
            lines.append(
                f"| {member} | {d['geometry']} | {d['mode']} | {preset} | {d['step']} "
                f"| {num(r, 'rmse_interior_m_yr'):.4f} | {num(r, 'mean_bias_interior_m_yr'):+.4f} "
                f"| {num(r, 'rmse_reach_interior_m_yr'):.4f} | {near_rmse:.4f} "
                f"| {num(r, f'be_rate_gis{first}_m_yr'):+.1f} | {num(r, f'be_rate_gis{last}_m_yr'):+.1f} |")
    lines += ["", f"## Rates on GIS {lo}-{hi}, m/yr (model LRR; target is the surveyed LOESS)", "",
              "| GIS | target | " + " | ".join(f"{m} {p}" for m in order for p in ("zeroBE", "edgeBE")
                                                if (m, p) in finals) + " |",
              "|---|---|" + "|".join("---" for m in order for p in ("zeroBE", "edgeBE")
                                     if (m, p) in finals) + "|"]
    for g in near:
        cells = [f"{finals[(m, p)]['lrr'].get(g, float('nan')):+.3f}"
                 for m in order for p in ("zeroBE", "edgeBE") if (m, p) in finals]
        lines.append(f"| {g} | {target.get(g, float('nan')):+.3f} | " + " | ".join(cells) + " |")
    out = EXPERIMENT_DIR / "RESULTS.md"
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines[:len(order) * 2 + 8]))
    print(f"\nwrote {out}")
    draw(finals, target)
    draw_compare(finals, target)
    return 0


def draw_compare(finals, target):
    """The figure to read first: the main approach (base geometry, edgeBE,
    compressed planform, the matrix run) against the solved n115 extension
    on GIS 1-90 only. One panel (Hannah, 2026-09-16: no difference or
    residual panels, no southern extension)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from site_layer.hat_figure_style import (C, C_1997, INK_MUTED, apply_style, caption,
                                  figsize, save, DOMAIN_AXIS_LABEL, town_bands)
    apply_style()
    base = finals.get(("base-asrun", "edgeBE"))
    base0 = finals.get(("base-asrun", "zeroBE"))
    ext = finals.get(("n115-asrun", "edgeBE"))
    if base is None or ext is None:
        print("comparison figure needs the base edgeBE matrix run and the solved n115 run")
        return
    gis = np.arange(1, 91)
    pick = lambda d: np.array([d["lrr"].get(g, np.nan) for g in gis])
    tgt = np.array([target.get(g, np.nan) for g in gis])
    b, e = pick(base), pick(ext)

    fig, ax = plt.subplots(figsize=figsize("double", aspect=0.5), constrained_layout=True)
    ax.axhline(0, color=INK_MUTED, lw=0.5)
    ax.plot(gis, tgt, color=C["REF"], lw=1.6, label="Observed (CoastSat), 1996–2010")
    if base0 is not None:
        ax.plot(gis, pick(base0), color=C["BASE"], lw=0.9, ls=":",
                label="Modelled 1996–2010, 90-domain reach, no boundary terms")
    ax.plot(gis, b, color=C["BASE"], lw=1.6,
            label="Modelled 1996–2010, 90-domain reach, boundary terms at domains 1 and 90")
    ax.plot(gis, e, color=C_1997, lw=1.3,
            label="Modelled 1996–2010, 115-domain reach (Pea Island extension), boundary terms at domains 1 and 115")
    ax.set_xlim(0.5, 90.5)
    ax.set_ylabel("shoreline change, m/yr")
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    town_bands(ax, where="top")
    fig.legend(*ax.get_legend_handles_labels(), loc="outside lower left", ncol=1)

    def rmse(y):
        r = (y - tgt)[1:-1]
        return float(np.sqrt(np.nanmean(r ** 2)))
    caption(fig, "Modelled 1996-2010 shoreline change rate (OLS slope through the annual "
                 "states, + seaward) on GIS 1-90 for the main approach (the edgeBE matrix "
                 "run: the 90-domain reach with extrapolated buffers and boundary source/sink "
                 "terms of +32.2 at GIS 1 and +10.0 m/yr at GIS 90) and for the Pea Island "
                 "extension (GIS 1-115 on the measured coast, its ends solved to +37.1 "
                 "and +41.7 m/yr), against the CoastSat target (LOESS-10, raw means on "
                 "GIS 1-10); the main approach without boundary terms is dotted. Compressed "
                 "planform, full_management, no groin. The two runs coincide from GIS 13 "
                 "to 78 and part only where the GIS 90 boundary term acted (GIS 80-90), where "
                 "the extension follows the run without boundary terms. Interior RMSE over GIS 2-89: "
                 f"main approach {rmse(b):.3f}, extension {rmse(e):.3f} m/yr.")
    fig_dir = EXPERIMENT_DIR / "figures"
    save(fig, fig_dir / "compare_main_vs_extension_gis1_90.png")
    print(f"wrote {fig_dir / 'compare_main_vs_extension_gis1_90.png'}")


def draw(finals, target):
    """The whole extended reach: the alongshore LRR of the solved runs
    against the target, a panel per offset mode, in house style."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer.hat_figure_style import (C, C_1997, INK_MUTED, apply_style, caption,
                                  figsize, save, _title, DOMAIN_AXIS_LABEL)
    apply_style()
    fig, axes = plt.subplots(2, 1, figsize=figsize("double", aspect=0.9), sharex=True,
                             constrained_layout=True)
    styles = {"base": (C["BASE"], "Modelled, 90-domain reach"),
              "n115": (C_1997, "Modelled, 115-domain reach (Pea Island extension)")}
    ext_target = _target(extended=True)
    outside = ext_target[ext_target.index > 90]
    for ax, mode in zip(axes, ("asrun", "detrended")):
        ax.axhline(0, color=INK_MUTED, lw=0.5)
        ax.plot(target.index, target.values, color=C["REF"], lw=1.4, label="Observed (CoastSat), 1996–2010, surveyed reach")
        ax.plot(outside.index, outside.values, color=C["REF"], lw=1.0, ls="--",
                label="Observed (CoastSat), 1996–2010, extension")
        for geometry, (colour, label) in styles.items():
            member = f"{geometry}-{mode}" if not (geometry == "base" and mode == "asrun") else "base-asrun"
            d = finals.get((member, "edgeBE")) or finals.get((member, "zeroBE"))
            if d is None:
                continue
            first, last = GEOMETRIES[d["geometry"]]
            terms = (f"boundary terms at domains {first} and {last}"
                     if d["preset"] == "edgeBE" else "no boundary terms")
            ax.plot(d["lrr"].index, d["lrr"].values, color=colour, lw=1.0,
                    label=f"{label}, {terms}")
        ax.axvspan(90.5, 116, color="0.94", zorder=0)
        ax.set_ylabel("shoreline change, m/yr")
        _title(ax, list(("asrun", "detrended")).index(mode),
               "compressed planform (as run)" if mode == "asrun" else "detrended planform")
    axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)
    axes[-1].set_xlim(0.5, 116)
    fig.legend(*axes[0].get_legend_handles_labels(), loc="outside lower left", ncol=1)
    caption(fig, "Modelled 1996-2010 shoreline change rate (OLS slope through the annual "
                 "states, + seaward) along the reach for the 90-domain reach and the "
                 "extended reach, each with its boundary source/sink terms solved against "
                 "the observed rate at its end domains, against the CoastSat rate (LOESS-10; "
                 "raw means on GIS 1-10): solid over the surveyed reach, dashed over the "
                 "extension, where the LOESS is one-sided at GIS 115 as it is at GIS 90 in "
                 "the 90-domain reach. Shaded: beyond GIS 90. (a) the compressed planform "
                 "every calibrated run uses; (b) the detrended planform at full strength, "
                 "note the axis.")
    fig_dir = EXPERIMENT_DIR / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    save(fig, fig_dir / "alongshore_rates_extension.png")
    print(f"wrote {fig_dir / 'alongshore_rates_extension.png'}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name, fn in (("stage0", cmd_stage0), ("check", cmd_check)):
        p = sub.add_parser(name)
        p.add_argument("--overwrite", action="store_true")
        p.add_argument("--dry-run", action="store_true")
        if name == "stage0":
            p.add_argument("--members", nargs="+", choices=sorted(MEMBERS), default=None,
                           help="only these members (default: every member but base-check)")
        p.set_defaults(fn=fn)
    p = sub.add_parser("probe")
    p.add_argument("--member", required=True, choices=sorted(MEMBERS))
    p.add_argument("--step", type=int, required=True)
    p.add_argument("--override", required=True, help='e.g. "1=32.2,115=12.0"')
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    p.set_defaults(fn=cmd_probe)
    p = sub.add_parser("next")
    p.add_argument("--member", required=True, choices=sorted(MEMBERS))
    p.set_defaults(fn=cmd_next)
    p = sub.add_parser("solve")
    p.add_argument("--member", required=True, choices=sorted(MEMBERS))
    p.add_argument("--max-steps", type=int, default=4)
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    p.set_defaults(fn=cmd_solve)
    p = sub.add_parser("score")
    p.set_defaults(fn=cmd_score)
    a = ap.parse_args(argv)
    return a.fn(a)


if __name__ == "__main__":
    sys.exit(main())
