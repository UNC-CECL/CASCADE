r"""
HAT_offset_scale_wave_tuning.py -- island-offset scale against wave-climate tuning
==============================================================================
THE QUESTION (Hannah, 2026-09-24). BRIE's shoreline is in metres, and so is
the island-offset file, but every calibrated run hands BRIE the offset divided
by ten (`offset_mode: asrun`, a units error). Put back at full scale, can the
model be re-tuned through its wave climate to match the /10 runs' skill?

TWO SWEEPS, ONE STUDY
    wave_height  Hs 0.5-3.0 m at the default wave angles, every offset
                 scale, dune-line and shoreline sources
    wave_angle   wave asymmetry 0.5-0.8 x high-angle fraction 0.1-0.4 at
                 Hs 1.0 m, every offset scale, dune-line source
    fixed        1996-2010, zeroBE, full_management, no groin, relocations
                 off, base geometry, Tp 8 s

OFFSET SCALE -- the folder label, and the value the runner reads
    div10              HAT_OFFSET_MODE=asrun      offset / 10 (every calibrated run)
    metres             HAT_OFFSET_MODE=metres     the measurement as is
    metres-detrended   HAT_OFFSET_MODE=detrended  metres, linear trend removed

LAYOUT (output/raw_runs/experiments/<STUDY>/)
    README.md                     question, design, how to read labels, results
    tables/                       written by `score`, every setting a column;
                                  observed_target.csv: the target's mean, sd
                                  and flat-line RMSE
    figures/<sweep or combined>/  written by HAT_plot_offset_scale_wave_tuning.py
    logs/<sweep>/<scale>_<source>/Hs<h>_asymmetry<a>_highangle<f>.log
    logs/drivers/                 this script's own console logs
    runs_<sweep>/<scale>_<source>/1996_2010/zeroBE/<run_name>/

    The run folder's tag is <STUDY>/runs_<sweep>/<scale>_<source>, the
    runner's three-level maximum. The run NAME is the runner's and leaves out
    whatever is at its default (no offset token for div10, no asym token at
    0.7, no ahf token at 0.1); the folder, the log name and the tables always
    carry every setting.

USAGE
    python HAT_offset_scale_wave_tuning.py run wave_height [--scales ...] [--hs ...]
    python HAT_offset_scale_wave_tuning.py run wave_angle
    python HAT_offset_scale_wave_tuning.py score
    (--jobs N, --dry-run, --overwrite on `run`)
==============================================================================
"""
from __future__ import annotations

import argparse
import ctypes
import json
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from itertools import product
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve()
PROJECT_ROOT = next(_p for _p in _HERE.parents
                    if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

HINDCAST = PROJECT_ROOT / "scripts" / "hatteras_ms" / "HAT_hindcast_1984_2024.py"
RAW_RUNS = PROJECT_ROOT / "output" / "raw_runs"
STUDY = "2026-09-24-island-offset-scale-wave-tuning"
STUDY_DIR = RAW_RUNS / "experiments" / STUDY
TABLES_DIR = STUDY_DIR / "tables"
LOGS_DIR = STUDY_DIR / "logs"
RUN_TIMEOUT_S = 3600

PERIOD = 1996
PRESET = "zeroBE"
SCENARIO = "full_management"

# folder label -> the runner's HAT_OFFSET_MODE value
SCALES = {"div10": "asrun", "metres": "metres", "metres-detrended": "detrended"}
MODE_TO_SCALE = {v: k for k, v in SCALES.items()}
SOURCES = ("duneline", "shoreline")

# The calibration defaults a run is at unless a sweep moves it.
DEFAULT_ASYMMETRY = 0.7
DEFAULT_HIGH_FRACTION = 0.1

SWEEPS = {
    "wave_height": dict(
        scales=tuple(SCALES), sources=SOURCES,
        hs=(1.0, 1.5, 2.0, 2.5, 3.0),
        asymmetry=(DEFAULT_ASYMMETRY,), high_fraction=(DEFAULT_HIGH_FRACTION,),
        # Run after the first grid, as their own commands (see README):
        #   --scales metres metres-detrended --hs 0.5 0.75
        #   --scales div10 --hs 0.5
    ),
    "wave_angle": dict(
        scales=tuple(SCALES), sources=("duneline",),
        hs=(1.0,),
        asymmetry=(0.5, 0.6, 0.7, 0.8), high_fraction=(0.1, 0.2, 0.3, 0.4),
        # metres was still improving at 0.4, so extended (2026-09-24):
        #   --scales metres --high-fraction 0.45 0.5
    ),
}


def num(x):
    """1.0 -> '1.0', 0.75 -> '0.75': one spelling for every label."""
    s = f"{float(x):g}"
    return s if "." in s else s + ".0"


def member(scale, source):
    return f"{scale}_{source}"


def tag(sweep, scale, source):
    return f"{STUDY}/runs_{sweep}/{member(scale, source)}"


def log_path(sweep, scale, source, hs, asym, ahf):
    return (LOGS_DIR / sweep / member(scale, source)
            / f"Hs{num(hs)}_asymmetry{num(asym)}_highangle{num(ahf)}.log")


def run_env(sweep, scale, source, hs, asym, ahf, overwrite=False):
    """The environment one run reads, built the way HAT_run_all builds it:
    every HAT_* variable named here, none inherited from the shell."""
    env = {k: v for k, v in os.environ.items() if not k.startswith("HAT_")}
    env.update({
        "HAT_IGNORE_SETTINGS": "1",
        "HAT_START_YEAR": str(PERIOD),
        "HAT_SOURCE_SINK_PRESET": PRESET,
        "HAT_SCENARIO": SCENARIO,
        "HAT_RELOCATIONS": "false",
        "HAT_GROIN_ENABLED": "false",
        "HAT_OFFSET_MODE": SCALES[scale],
        "HAT_ISLAND_OFFSET_SOURCE": source,
        "HAT_HS": f"{hs}",
        "HAT_WAVE_ASYMMETRY": f"{asym}",
        "HAT_WAVE_ANGLE_HIGH_FRACTION": f"{ahf}",
        "HAT_RUN_KIND": "experiment",
        "HAT_RUN_TAG": tag(sweep, scale, source),
        "HAT_OVERWRITE": "true" if overwrite else "false",
        "HAT_SAVE_MODEL_STATE": "false",
        "HAT_MAKE_GIFS": "false",
        "MPLBACKEND": "Agg",
        "PYTHONIOENCODING": "utf-8",
    })
    return env


def launch(sweep, scale, source, hs, asym, ahf, overwrite=False, dry_run=False):
    label = (f"{sweep} {member(scale, source)} Hs {num(hs)} "
             f"asymmetry {num(asym)} high-angle {num(ahf)}")
    log = log_path(sweep, scale, source, hs, asym, ahf)
    # A clean log means the run exists; the runner would refuse it anyway,
    # and that refusal used to read as a failure.
    if not overwrite and log.is_file() and "Traceback" not in log.read_text(
            encoding="utf-8", errors="replace"):
        print(f"skip {label}: already run", flush=True)
        return True
    if dry_run:
        print(f"would run: {label}")
        return True
    log.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    proc = subprocess.run([sys.executable, str(HINDCAST)],
                          env=run_env(sweep, scale, source, hs, asym, ahf, overwrite),
                          cwd=str(PROJECT_ROOT), capture_output=True, text=True,
                          encoding="utf-8", errors="replace", timeout=RUN_TIMEOUT_S)
    log.write_text((proc.stdout or "") + "\n--- STDERR ---\n" + (proc.stderr or ""),
                   encoding="utf-8")
    minutes = (time.perf_counter() - t0) / 60
    if proc.returncode != 0:
        print(f"FAILED {label} after {minutes:.1f} min (exit {proc.returncode}): "
              f"{stop_reason(log)}; log {log.relative_to(STUDY_DIR)}", flush=True)
        return False
    print(f"done {label} in {minutes:.1f} min", flush=True)
    return True


def stop_reason(log):
    """Why a run has no score, read from its log: the drowning, or the error."""
    text = log.read_text(encoding="utf-8", errors="replace")
    m = re.search(r"Model stopped at year (\d+)", text)
    if m:
        kinds = sorted(set(re.findall(r"Barrier has (\w+) DROWNED", text)))
        return f"barrier drowned ({'/'.join(k.lower() for k in kinds)}) in year {m.group(1)}"
    lines = [ln for ln in text.splitlines() if ln.strip()]
    return lines[-1][:200] if lines else "empty log"


def keep_awake():
    """A sleeping machine kills a sweep mid-run; hold it awake while we live."""
    if os.name == "nt":
        ctypes.windll.kernel32.SetThreadExecutionState(0x80000000 | 0x00000001)


def cmd_run(a):
    spec = SWEEPS[a.sweep]
    cells = list(product([a.sweep], a.scales or spec["scales"], spec["sources"],
                         a.hs or spec["hs"], a.asymmetry or spec["asymmetry"],
                         a.high_fraction or spec["high_fraction"]))
    keep_awake()
    print(f"{a.sweep}: {len(cells)} runs, {a.jobs} at a time", flush=True)
    with ThreadPoolExecutor(max_workers=a.jobs) as pool:
        ok = list(pool.map(lambda c: launch(*c, overwrite=a.overwrite,
                                            dry_run=a.dry_run), cells))
    print(f"{sum(ok)} of {len(cells)} succeeded")
    return 0 if all(ok) else 1


# =============================================================================
# SCORE
# =============================================================================

# THE TARGET AND THE ALONGSHORE SCORES
#
# The runner scores RMSE and bias; these add how much of the observed
# alongshore variation a run explains. The target is rebuilt exactly as the
# runner builds it (section 8 of HAT_hindcast_1984_2024.py), and every run's
# RMSE is recomputed from it and checked against the runner's, so the two
# sets of scores are provably on the same target.

def coastsat_target(start=PERIOD):
    """The CoastSat LRR target, GIS 1-90, as the runner builds it for a start
    year (section 8 of the runner: LOESS at 10 domains, the southern 10 raw).
    Shared with scripts/sensitivity_analysis/natural_wave_sensitivity.py."""
    from site_layer.hat_observed_rates import COASTSAT_LRR_ROOT
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS, HATTERAS_PERIODS
    from cascade_pipeline.hindcast import build_target_table
    from cascade_pipeline.coastsat_loess import (CoastSatDataset, LoessConfig,
                                                 build_coastsat_series)
    window = f"{start}_{HATTERAS_PERIODS[start]['end_year']}"
    ds = CoastSatDataset(label=f"CoastSat LRR ({window.replace('_', '-')})",
                         period_start=start,
                         csv_path=str(COASTSAT_LRR_ROOT / window / "transect_lrr_full.csv"))
    cfg = LoessConfig(window_domains=(10,), skip_southern_domains=10)
    cs = build_coastsat_series([ds], active_period_start=start, loess_config=cfg,
                               domains=HATTERAS_DOMAINS)[0]
    return build_target_table(cs, cfg, HATTERAS_DOMAINS, 10).set_index(
        "gis_domain")["target_lrr_m_yr"]


def interior(series):
    from site_layer.hatteras_site_config import SCORE_INTERIOR_GIS
    lo, hi = SCORE_INTERIOR_GIS
    return series.loc[lo:hi]


def run_rates(run_dir):
    import pandas as pd
    return pd.read_csv(Path(run_dir) / "tables" / "shoreline_change_rate.csv"
                       ).set_index("gis_domain")["lrr_m_yr"]


def alongshore_scores(rates, target):
    """How much of the observed alongshore variation a run explains.

    variance_explained          1 - sum((m - o)^2) / sum((o - mean o)^2).
                                1 is perfect; 0 is no better than a flat
                                line at the observed mean; negative is worse.
                                Bias counts against it.
    pattern_variance_explained  the same with each series' own mean removed
                                first: the pattern alone, bias forgiven.
    r_alongshore                correlation of the two alongshore series
    sd_ratio                    model sd / observed sd: below 1 the model
                                varies less along the island than observed
    """
    m, o = interior(rates), interior(target)
    sst = float(((o - o.mean()) ** 2).sum())
    return {
        "variance_explained": 1 - float(((m - o) ** 2).sum()) / sst,
        "pattern_variance_explained":
            1 - float((((m - m.mean()) - (o - o.mean())) ** 2).sum()) / sst,
        "r_alongshore": float(np.corrcoef(m, o)[0, 1]),
        "model_sd_m_yr": float(m.std(ddof=0)),
        "sd_ratio": float(m.std(ddof=0) / o.std(ddof=0)),
        "_rmse": float(np.sqrt(((m - o) ** 2).mean())),
    }


def _settings_from_log_name(path):
    m = re.fullmatch(r"Hs([\d.]+)_asymmetry([\d.]+)_highangle([\d.]+)\.log", path.name)
    return tuple(float(g) for g in m.groups())


def cmd_score(a):
    import pandas as pd
    from cascade_pipeline.run_registry import load_run_index, rebuild_run_index

    rebuild_run_index(RAW_RUNS)
    target = coastsat_target()
    index = load_run_index(RAW_RUNS / "run_index.csv")
    index = index[index["tag"].astype(str).str.startswith(STUDY + "/")
                  & (index["status"] == "current")]
    by_dir = {}
    for _, r in index.iterrows():
        run_dir = (RAW_RUNS / "experiments" / r["tag"]
                   / f"{r['start_year']}_{r['end_year']}" / r["source_sink_preset"]
                   / r["run_name"])
        by_dir[run_dir] = r

    records = []
    for log in sorted(LOGS_DIR.glob("wave_*/*/*.log")):
        sweep, mem = log.parts[-3], log.parts[-2]
        scale, source = mem.rsplit("_", 1)
        hs, asym, ahf = _settings_from_log_name(log)
        rec = {
            "sweep": sweep,
            "offset_scale": scale,
            "offset_mode_setting": SCALES[scale],
            "offset_source": source,
            "Hs_m": hs, "wave_asymmetry": asym, "wave_angle_high_fraction": ahf,
            "is_control": scale == "div10",
            "period": "1996-2010", "source_sink_preset": PRESET,
            "log": str(log.relative_to(STUDY_DIR)).replace("\\", "/"),
        }
        # the run this log belongs to: the one whose metadata says these settings
        match = None
        for run_dir, r in by_dir.items():
            if run_dir.parent.parent.parent.name != mem or run_dir.parts[-5] != f"runs_{sweep}":
                continue
            md = json.loads((run_dir / f"{r['run_name']}_run_metadata.json")
                            .read_text(encoding="utf-8"))
            w = md["wave climate"]
            if (float(w["wave_height_m"]), float(w["wave_asymmetry"]),
                    float(w["wave_angle_high_frac"])) == (hs, asym, ahf):
                md_mode = md["scenario"]["shoreline offset"]
                md_version = md["identity"]["island_offset_version"]
                if md_mode != SCALES[scale] or not md_version.startswith(source + "/"):
                    raise ValueError(f"{run_dir}: metadata says {md_mode!r} / "
                                     f"{md_version!r}, folder says {mem!r}")
                match = (run_dir, r, md)
                break
        if match:
            run_dir, r, md = match
            scores = alongshore_scores(run_rates(run_dir), target)
            if not np.isclose(scores.pop("_rmse"), float(r["rmse_interior_m_yr"]), rtol=1e-3):
                raise ValueError(f"{run_dir}: RMSE against the rebuilt target does "
                                 f"not match the runner's -- not the same target")
            rec.update(scores)
            rec.update({
                "status": "scored",
                "island_offset_version": md["identity"]["island_offset_version"],
                "wave_period_s": float(md["wave climate"]["wave_period_s"]),
                "mean_bias_interior_m_yr": r["mean_bias_interior_m_yr"],
                "rmse_interior_m_yr": r["rmse_interior_m_yr"],
                "endpoint_mean_bias_interior_m_yr": r["endpoint_mean_bias_interior_m_yr"],
                "endpoint_rmse_interior_m_yr": r["endpoint_rmse_interior_m_yr"],
                "lrr_r2_median": r["lrr_r2_median"],
                "roads_drowned": r["roads_drowned"],
                "run_name": r["run_name"],
                "run_dir": str(run_dir.relative_to(STUDY_DIR)).replace("\\", "/"),
            })
        else:
            rec["status"] = stop_reason(log)
        records.append(rec)

    cols = ["sweep", "offset_scale", "offset_mode_setting", "offset_source",
            "island_offset_version", "is_control", "Hs_m", "wave_period_s",
            "wave_asymmetry", "wave_angle_high_fraction", "period",
            "source_sink_preset", "status", "mean_bias_interior_m_yr",
            "rmse_interior_m_yr", "variance_explained",
            "pattern_variance_explained", "r_alongshore", "model_sd_m_yr",
            "sd_ratio", "endpoint_mean_bias_interior_m_yr",
            "endpoint_rmse_interior_m_yr", "lrr_r2_median", "roads_drowned",
            "run_name", "run_dir", "log"]
    out = (pd.DataFrame(records).reindex(columns=cols)
           .sort_values(["sweep", "offset_scale", "offset_source", "Hs_m",
                         "wave_asymmetry", "wave_angle_high_fraction"]))
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    out.to_csv(TABLES_DIR / "all_runs.csv", index=False)
    for sweep in SWEEPS:
        out[out.sweep == sweep].to_csv(TABLES_DIR / f"{sweep}_sweep.csv", index=False)
    t = interior(target)
    pd.DataFrame([{
        "target": "CoastSat LRR 1996-2010, LOESS 10 domains",
        "domains": "GIS 2-89", "n_domains": len(t),
        "mean_m_yr": t.mean(), "sd_m_yr": t.std(ddof=0),
        "flat_line_rmse_m_yr": t.std(ddof=0),
        "note": "flat line = the observed mean at every domain; "
                "variance_explained = 1 - (rmse / flat_line_rmse)^2",
    }]).to_csv(TABLES_DIR / "observed_target.csv", index=False)
    print(out.groupby(["sweep", "status"]).size().to_string())
    return 0


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("sweep", choices=tuple(SWEEPS))
    r.add_argument("--scales", nargs="+", choices=tuple(SCALES))
    r.add_argument("--hs", nargs="+", type=float)
    r.add_argument("--asymmetry", nargs="+", type=float)
    r.add_argument("--high-fraction", nargs="+", type=float)
    r.add_argument("--jobs", type=int, default=6)
    r.add_argument("--dry-run", action="store_true")
    r.add_argument("--overwrite", action="store_true")
    r.set_defaults(func=cmd_run)
    s = sub.add_parser("score")
    s.set_defaults(func=cmd_score)
    a = p.parse_args()
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
