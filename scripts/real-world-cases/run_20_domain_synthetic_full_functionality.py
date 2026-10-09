#!/usr/bin/env python3
"""Run and visualize a 20-domain synthetic RoadwayManager integration test.

This is a full CASCADE/Barrier3D run, not a prescribed trajectory. All domains
use the same established synthetic Barrier3D elevation, dune, growth, and storm
inputs. The only externally requested actions are the documented test-matrix
events below; all storm, overwash, dune, shoreline, roadway, and drowning
responses are calculated by the model.

The script deliberately writes to a new output directory and refuses to replace
an existing run. It does not edit RoadwayManager or any other model source.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import sys
from pathlib import Path

import matplotlib
import numpy as np
import yaml
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


SOURCE_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = SOURCE_ROOT.parent
sys.path.insert(0, str(SOURCE_ROOT))

from cascade import Cascade  # noqa: E402


DOMAIN_COUNT = 20
TIME_STEP_COUNT = 200
SETBACK_PATTERN_M = [10.0, 20.0, 21.0, 30.0]
ROAD_SETBACKS_M = SETBACK_PATTERN_M * 5
ROAD_SETBACK_TRIGGER_M = 20.0
NOURISHMENT_VOLUME_M3_PER_M = 300.0
AUTOMATIC_NOURISHMENT_INTERVAL_YR = 20
MANUAL_NOURISHMENT_YEAR = 1
HISTORICAL_RELOCATION_YEAR = 1
HISTORICAL_RELOCATION_SETBACK_M = 30.0

GROUPS = (
    ("Control", range(0, 4)),
    ("Manual nourishment", range(4, 8)),
    ("Automatic nourishment", range(8, 12)),
    ("Historical relocation", range(12, 16)),
    ("Combined", range(16, 20)),
)
GROUP_FOR_DOMAIN = {
    domain: group for group, domains in GROUPS for domain in domains
}
GROUP_COLORS = {
    "Control": "#636363",
    "Manual nourishment": "#2171b5",
    "Automatic nourishment": "#6a51a3",
    "Historical relocation": "#d95f0e",
    "Combined": "#238b45",
}

INPUT_ROOT = SOURCE_ROOT / "tests" / "test_human_dynamics"
INPUT_FILES = (
    "StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
    "b3d_pt45_8750yrs_low-elevations.csv",
    "pathways-dunes.npy",
    "growthparam_1000dam.npy",
    "nourishment-parameters.yaml",
)
DEFAULT_OUTPUT_DIR = (
    WORKSPACE_ROOT
    / "CASCADE"
    / "output"
    / "roadway_functionalities_full_synthetic"
    / "20_domains_200_years"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Fresh directory for this run (existing directories are rejected).",
    )
    parser.add_argument(
        "--skip-gif",
        action="store_true",
        help="Run the model and make static outputs without rendering the GIF.",
    )
    return parser.parse_args()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def prepare_fresh_run(output_dir: Path) -> Path:
    if output_dir.exists():
        raise FileExistsError(
            f"Refusing to replace an existing run directory: {output_dir}"
        )
    input_dir = output_dir / "inputs"
    input_dir.mkdir(parents=True)
    for filename in INPUT_FILES:
        source = INPUT_ROOT / filename
        if not source.is_file():
            raise FileNotFoundError(f"Required synthetic input is missing: {source}")
        shutil.copy2(source, input_dir / filename)
    return input_dir


def validate_inputs(input_dir: Path) -> dict:
    storms = np.load(
        input_dir / "StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy"
    )
    dunes = np.load(input_dir / "pathways-dunes.npy")
    growth = np.load(input_dir / "growthparam_1000dam.npy")
    elevation = np.loadtxt(
        input_dir / "b3d_pt45_8750yrs_low-elevations.csv", delimiter=","
    )
    with (input_dir / "nourishment-parameters.yaml").open() as stream:
        parameters = yaml.safe_load(stream)

    storm_start = int(parameters["StormStart"])
    if storm_start >= TIME_STEP_COUNT:
        raise ValueError(
            f"StormStart={storm_start} disables storms in this {TIME_STEP_COUNT}-year run."
        )
    if int(np.nanmax(storms[:, 0])) < TIME_STEP_COUNT - 1:
        raise ValueError("The storm series does not cover the requested run duration.")
    if dunes.size < int(parameters["BarrierLength"]):
        raise ValueError("The dune input is shorter than the synthetic domain.")
    if growth.size < int(parameters["BarrierLength"]):
        raise ValueError("The growth-parameter input is shorter than the domain.")
    if len(ROAD_SETBACKS_M) != DOMAIN_COUNT:
        raise ValueError("The road-setback matrix must contain exactly 20 values.")

    return {
        "storm_shape": list(storms.shape),
        "storm_start": storm_start,
        "storm_last_year": int(np.nanmax(storms[:, 0])),
        "elevation_shape": list(elevation.shape),
        "dune_shape": list(dunes.shape),
        "growth_shape": list(growth.shape),
    }


def domain_matrix() -> list[dict]:
    rows = []
    for domain in range(DOMAIN_COUNT):
        group = GROUP_FOR_DOMAIN[domain]
        rows.append(
            {
                "domain": domain,
                "group": group,
                "initial_road_setback_m": ROAD_SETBACKS_M[domain],
                "road_setback_trigger_m": ROAD_SETBACK_TRIGGER_M,
                "manual_nourishment_year": (
                    MANUAL_NOURISHMENT_YEAR
                    if group in {"Manual nourishment", "Combined"}
                    else None
                ),
                "automatic_nourishment_interval_yr": (
                    AUTOMATIC_NOURISHMENT_INTERVAL_YR
                    if group in {"Automatic nourishment", "Combined"}
                    else None
                ),
                "historical_relocation_year": (
                    HISTORICAL_RELOCATION_YEAR
                    if group in {"Historical relocation", "Combined"}
                    else None
                ),
                "historical_relocation_target_setback_m": (
                    HISTORICAL_RELOCATION_SETBACK_M
                    if group in {"Historical relocation", "Combined"}
                    else None
                ),
            }
        )
    return rows


def cascade_parameters() -> dict:
    intervals = []
    for domain in range(DOMAIN_COUNT):
        if GROUP_FOR_DOMAIN[domain] in {"Automatic nourishment", "Combined"}:
            intervals.append(AUTOMATIC_NOURISHMENT_INTERVAL_YR)
        else:
            intervals.append(None)

    return {
        "storm_file": "StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
        "elevation_file": "b3d_pt45_8750yrs_low-elevations.csv",
        "dune_file": "pathways-dunes.npy",
        "parameter_file": "nourishment-parameters.yaml",
        "alongshore_section_count": DOMAIN_COUNT,
        "time_step_count": TIME_STEP_COUNT,
        "num_cores": 1,
        "roadway_management_module": [True] * DOMAIN_COUNT,
        "beach_nourishment_module": [False] * DOMAIN_COUNT,
        "community_economics_module": False,
        "outwash_module": [False] * DOMAIN_COUNT,
        "alongshore_transport_module": True,
        "road_setback": ROAD_SETBACKS_M,
        "road_setback_trigger": [ROAD_SETBACK_TRIGGER_M] * DOMAIN_COUNT,
        "nourishment_interval": intervals,
        "nourishment_volume": [NOURISHMENT_VOLUME_M3_PER_M] * DOMAIN_COUNT,
    }


def write_manifest(
    output_dir: Path,
    input_dir: Path,
    input_validation: dict,
    status: str,
    completed_year: int | None = None,
    stop_reason: str | None = None,
) -> None:
    supplied = cascade_parameters()
    manifest = {
        "status": status,
        "run_name": "roadway_functionalities_20_domain_synthetic",
        "source_root": str(SOURCE_ROOT),
        "model_file": str(SOURCE_ROOT / "cascade" / "roadway_manager.py"),
        "requested_time_steps": TIME_STEP_COUNT,
        "completed_year": completed_year,
        "stop_reason": stop_reason,
        "domain_count": DOMAIN_COUNT,
        "domain_matrix": domain_matrix(),
        "cascade_arguments_changed_from_defaults": supplied,
        "defaults_left_unchanged": {
            "wave_height_m": 1,
            "wave_period_s": 7,
            "wave_asymmetry": 0.8,
            "wave_angle_high_fraction": 0.2,
            "bay_depth_m": 3.0,
            "background_slope": 0.001,
            "berm_elevation_m_NAVD88": 1.9,
            "MHW_m_NAVD88": 0.46,
            "beach_slope": 0.04,
            "RSLR_m_per_year": 0.004,
            "RSLR_constant": True,
            "background_erosion_m_per_year": 0.0,
            "min_dune_growth_rate": 0.25,
            "max_dune_growth_rate": 0.65,
            "road_elevation_m_MHW": 1.7,
            "road_width_m": 30,
            "dune_design_elevation_m_MHW": 3.7,
            "dune_minimum_elevation_m_MHW": 2.2,
        },
        "important_module_choices": {
            "RoadwayManager": "on for all 20 domains",
            "BeachDuneManager": "off for all 20 domains",
            "BRIE_alongshore_transport": "on (CASCADE default)",
            "Outwasher": "off to keep this a RoadwayManager integration test",
            "CHOM": "off (CASCADE default)",
        },
        "test_event_convention": (
            "Manual nourishment and historical relocation are queued before the "
            "first Cascade.update and are recorded at output year/index 1."
        ),
        "input_validation": input_validation,
        "input_sha256": {
            filename: file_sha256(input_dir / filename) for filename in INPUT_FILES
        },
        "output_directory": str(output_dir),
    }
    with (output_dir / "manifest.json").open("w") as stream:
        json.dump(manifest, stream, indent=2)


def run_model(input_dir: Path) -> tuple[Cascade, int, str]:
    cascade = Cascade(
        str(input_dir),
        name="roadway_functionalities_20_domain_synthetic",
        **cascade_parameters(),
    )

    completed_year = 0
    stop_reason = "completed requested duration"
    for year in range(1, TIME_STEP_COUNT):
        if year == MANUAL_NOURISHMENT_YEAR:
            requests = [0] * DOMAIN_COUNT
            for domain in (*range(4, 8), *range(16, 20)):
                requests[domain] = 1
            cascade.nourish_now = requests

        if year == HISTORICAL_RELOCATION_YEAR:
            for domain in range(12, 20):
                cascade.roadways[domain].request_historical_relocation(
                    HISTORICAL_RELOCATION_SETBACK_M
                )

        cascade.update()
        available = min(len(barrier.x_s_TS) for barrier in cascade.barrier3d) - 1
        completed_year = max(completed_year, available)
        if cascade.b3d_break:
            drowned = [
                domain
                for domain, barrier in enumerate(cascade.barrier3d)
                if barrier.drown_break
            ]
            stop_reason = f"Barrier3D drowning; domains={drowned}"
            print(f"Stopped at year {completed_year}: {stop_reason}", flush=True)
            break
        if year == 1 or year % 10 == 0 or year == TIME_STEP_COUNT - 1:
            active_roads = sum(not bool(value) for value in cascade.road_break)
            print(
                f"Completed year {year:3d}/{TIME_STEP_COUNT - 1}; "
                f"active roads={active_roads}/{DOMAIN_COUNT}",
                flush=True,
            )

    return cascade, completed_year, stop_reason


def pad(values, count: int, fill=np.nan, dtype=float) -> np.ndarray:
    result = np.full(count, fill, dtype=dtype)
    array = np.asarray(values)
    copied = min(count, array.size)
    if copied:
        result[:copied] = array[:copied]
    return result


def extract_results(cascade: Cascade, completed_year: int) -> dict[str, np.ndarray]:
    count = completed_year + 1
    shape = (DOMAIN_COUNT, count)
    arrays = {
        "shoreline_m": np.full(shape, np.nan),
        "shoreline_change_cells": np.full(shape, np.nan),
        "dune_grid_m": np.full(shape, np.nan),
        "dune_crest_mean_m_MHW": np.full(shape, np.nan),
        "mean_interior_height_m_MHW": np.full(shape, np.nan),
        "beach_width_m": np.full(shape, np.nan),
        "road_setback_m": np.full(shape, np.nan),
        "road_elevation_m_MHW": np.full(shape, np.nan),
        "storm_count": np.zeros(shape, dtype=int),
        "overwash_flux_m3_per_m": np.full(shape, np.nan),
        "nourishment": np.zeros(shape, dtype=bool),
        "nourishment_volume_m3_per_m": np.zeros(shape),
        "dune_migration_on": np.zeros(shape, dtype=bool),
        "dunes_rebuilt": np.zeros(shape, dtype=bool),
        "rebuild_disabled_by_setback": np.zeros(shape, dtype=bool),
        "triggered_relocation": np.zeros(shape, dtype=bool),
        "historical_relocation_requested": np.zeros(shape, dtype=bool),
        "forced_relocation": np.zeros(shape, dtype=bool),
        "incomplete_relocation": np.zeros(shape, dtype=bool),
    }

    for domain, (barrier, manager) in enumerate(
        zip(cascade.barrier3d, cascade.roadways)
    ):
        shoreline = pad(np.asarray(barrier.x_s_TS) * 10.0, count)
        shoreline_change = pad(barrier.ShorelineChangeTS, count)
        beach = pad(manager.beach_width, count)
        initial_dune_line = shoreline[0] + beach[0]
        dune_grid = initial_dune_line + np.cumsum(
            np.nan_to_num(-shoreline_change, nan=0.0)
        ) * 10.0

        dune_domain = np.asarray(barrier.DuneDomain[:count], dtype=float)
        dune_crest = (
            np.mean(np.max(dune_domain, axis=2), axis=1) * 10.0
            + barrier.BermEl * 10.0
        )
        n_dune = min(count, dune_crest.size)

        road_width = pad(manager._road_width_TS, count)
        road_setback = pad(manager._road_setback_TS, count)
        road_elevation = pad(manager._road_ele_TS, count)
        road_missing = road_width <= 0
        road_setback[road_missing] = np.nan
        road_elevation[road_missing] = np.nan

        arrays["shoreline_m"][domain] = shoreline
        arrays["shoreline_change_cells"][domain] = shoreline_change
        arrays["dune_grid_m"][domain] = dune_grid
        arrays["dune_crest_mean_m_MHW"][domain, :n_dune] = dune_crest[:n_dune]
        arrays["mean_interior_height_m_MHW"][domain] = pad(
            np.asarray(barrier.h_b_TS) * 10.0, count
        )
        arrays["beach_width_m"][domain] = beach
        arrays["road_setback_m"][domain] = road_setback
        arrays["road_elevation_m_MHW"][domain] = road_elevation
        arrays["storm_count"][domain] = pad(
            barrier._StormCount, count, fill=0, dtype=int
        )
        arrays["overwash_flux_m3_per_m"][domain] = pad(barrier.QowTS, count)
        arrays["nourishment"][domain] = pad(
            manager.nourishment_TS, count, fill=False, dtype=bool
        )
        arrays["nourishment_volume_m3_per_m"][domain] = pad(
            manager.nourishment_volume_TS, count
        )
        migration = pad(manager.dune_migration_on, count)
        arrays["dune_migration_on"][domain] = np.nan_to_num(
            migration, nan=1.0
        ).astype(bool)
        arrays["dunes_rebuilt"][domain] = pad(
            manager._dunes_rebuilt_TS, count, fill=False, dtype=bool
        )
        arrays["rebuild_disabled_by_setback"][domain] = pad(
            manager.road_dune_rebuild_disabled_TS,
            count,
            fill=False,
            dtype=bool,
        )
        arrays["triggered_relocation"][domain] = pad(
            manager.triggered_relocation_TS, count, fill=False, dtype=bool
        )
        arrays["historical_relocation_requested"][domain] = pad(
            manager.historical_relocation_requested_TS,
            count,
            fill=False,
            dtype=bool,
        )
        arrays["forced_relocation"][domain] = pad(
            manager.forced_relocation_TS, count, fill=False, dtype=bool
        )
        arrays["incomplete_relocation"][domain] = pad(
            manager.relocation_incomplete_TS, count, fill=False, dtype=bool
        )

    arrays["year"] = np.arange(count)
    arrays["domain"] = np.arange(DOMAIN_COUNT)
    arrays["group"] = np.asarray(
        [GROUP_FOR_DOMAIN[domain] for domain in range(DOMAIN_COUNT)]
    )
    return arrays


def years_where(values: np.ndarray) -> str:
    return ";".join(str(year) for year in np.flatnonzero(values))


def write_summary_csv(
    output_dir: Path, arrays: dict[str, np.ndarray], cascade: Cascade
) -> None:
    fields = [
        "domain",
        "group",
        "initial_setback_m",
        "final_recorded_setback_m",
        "road_active_at_end",
        "nourishment_years",
        "triggered_relocation_years",
        "historical_request_years",
        "forced_relocation_years",
        "incomplete_relocation_years",
        "dune_rebuild_years",
        "years_rebuild_disabled_by_setback",
        "total_storms",
        "cumulative_overwash_flux_m3_per_m",
        "final_shoreline_change_m",
        "final_mean_dune_crest_m_MHW",
        "final_mean_interior_height_m_MHW",
    ]
    with (output_dir / "domain_summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for domain in range(DOMAIN_COUNT):
            setbacks = arrays["road_setback_m"][domain]
            finite = setbacks[np.isfinite(setbacks)]
            shoreline = arrays["shoreline_m"][domain]
            writer.writerow(
                {
                    "domain": domain,
                    "group": GROUP_FOR_DOMAIN[domain],
                    "initial_setback_m": ROAD_SETBACKS_M[domain],
                    "final_recorded_setback_m": finite[-1] if finite.size else "",
                    "road_active_at_end": not bool(cascade.road_break[domain]),
                    "nourishment_years": years_where(arrays["nourishment"][domain]),
                    "triggered_relocation_years": years_where(
                        arrays["triggered_relocation"][domain]
                    ),
                    "historical_request_years": years_where(
                        arrays["historical_relocation_requested"][domain]
                    ),
                    "forced_relocation_years": years_where(
                        arrays["forced_relocation"][domain]
                    ),
                    "incomplete_relocation_years": years_where(
                        arrays["incomplete_relocation"][domain]
                    ),
                    "dune_rebuild_years": years_where(
                        arrays["dunes_rebuilt"][domain]
                    ),
                    "years_rebuild_disabled_by_setback": years_where(
                        arrays["rebuild_disabled_by_setback"][domain]
                    ),
                    "total_storms": int(np.nansum(arrays["storm_count"][domain])),
                    "cumulative_overwash_flux_m3_per_m": float(
                        np.nansum(arrays["overwash_flux_m3_per_m"][domain])
                    ),
                    "final_shoreline_change_m": shoreline[-1] - shoreline[0],
                    "final_mean_dune_crest_m_MHW": arrays[
                        "dune_crest_mean_m_MHW"
                    ][domain, -1],
                    "final_mean_interior_height_m_MHW": arrays[
                        "mean_interior_height_m_MHW"
                    ][domain, -1],
                }
            )


def write_validation_report(
    output_dir: Path,
    arrays: dict[str, np.ndarray],
    cascade: Cascade,
) -> list[dict]:
    """Evaluate expected test-matrix behavior without changing model state."""

    checks = []

    def record(name: str, passed: bool, observed: str, expected: str) -> None:
        checks.append(
            {
                "check": name,
                "passed": bool(passed),
                "observed": observed,
                "expected": expected,
            }
        )

    def event_years(key: str, domain: int) -> list[int]:
        return np.flatnonzero(arrays[key][domain]).tolist()

    record(
        "requested duration completed",
        int(arrays["year"][-1]) == TIME_STEP_COUNT - 1,
        str(int(arrays["year"][-1])),
        str(TIME_STEP_COUNT - 1),
    )
    total_storms = int(np.sum(arrays["storm_count"]))
    record(
        "active storms were simulated",
        total_storms > 0,
        str(total_storms),
        "> 0",
    )
    active_roads = sum(not bool(value) for value in cascade.road_break)
    record(
        "all roads remained active through final year",
        active_roads == DOMAIN_COUNT,
        str(active_roads),
        str(DOMAIN_COUNT),
    )

    control_events = sum(
        len(event_years("nourishment", domain)) for domain in range(0, 4)
    )
    record(
        "controls received no nourishment",
        control_events == 0,
        str(control_events),
        "0",
    )
    manual_observed = {
        domain: event_years("nourishment", domain) for domain in range(4, 8)
    }
    record(
        "manual nourishment occurred once at year 1",
        all(years == [1] for years in manual_observed.values()),
        json.dumps(manual_observed),
        "each domain [1]",
    )
    expected_automatic = list(range(20, TIME_STEP_COUNT, 20))
    automatic_observed = {
        domain: event_years("nourishment", domain) for domain in range(8, 12)
    }
    record(
        "automatic nourishment followed the 20-year interval",
        all(years == expected_automatic for years in automatic_observed.values()),
        json.dumps(automatic_observed),
        f"each domain {expected_automatic}",
    )
    expected_combined = [1] + list(range(21, TIME_STEP_COUNT, 20))
    combined_observed = {
        domain: event_years("nourishment", domain) for domain in range(16, 20)
    }
    record(
        "combined domains reset interval after manual nourishment",
        all(years == expected_combined for years in combined_observed.values()),
        json.dumps(combined_observed),
        f"each domain {expected_combined}",
    )

    historical_observed = {
        domain: event_years("historical_relocation_requested", domain)
        for domain in range(12, 20)
    }
    forced_observed = {
        domain: event_years("forced_relocation", domain)
        for domain in range(12, 20)
    }
    record(
        "historical relocation requests were recorded",
        all(years == [1] for years in historical_observed.values()),
        json.dumps(historical_observed),
        "each domain [1]",
    )
    record(
        "historical relocation requests completed successfully",
        all(years == [1] for years in forced_observed.values()),
        json.dumps(forced_observed),
        "each domain [1]",
    )
    incomplete_count = int(np.sum(arrays["incomplete_relocation"]))
    record(
        "no relocation was incomplete",
        incomplete_count == 0,
        str(incomplete_count),
        "0",
    )

    below_or_equal_domains = [0, 1, 4, 5, 8, 9]
    above_domains = [2, 3, 6, 7, 10, 11]
    below_states = arrays["rebuild_disabled_by_setback"][below_or_equal_domains, 1]
    above_states = arrays["rebuild_disabled_by_setback"][above_domains, 1]
    record(
        "dune rebuilding remained eligible at setbacks <= 20 m",
        not np.any(below_states),
        np.array2string(below_states.astype(int)),
        "all 0 (not disabled)",
    )
    record(
        "dune rebuilding was disabled at setbacks > 20 m",
        np.all(above_states),
        np.array2string(above_states.astype(int)),
        "all 1 (disabled)",
    )

    nourishment_domains = list(range(4, 12)) + list(range(16, 20))
    nourishment_response = []
    for domain in nourishment_domains:
        for year in event_years("nourishment", domain):
            nourishment_response.append(
                arrays["shoreline_m"][domain, year]
                < arrays["shoreline_m"][domain, year - 1]
            )
    record(
        "every nourishment event produced net shoreline progradation that year",
        bool(nourishment_response) and all(nourishment_response),
        f"{sum(nourishment_response)}/{len(nourishment_response)} events",
        f"{len(nourishment_response)}/{len(nourishment_response)} events",
    )

    triggered_count = int(np.sum(arrays["triggered_relocation"]))
    record(
        "native triggered relocation was observed",
        triggered_count > 0,
        str(triggered_count),
        "> 0 (diagnostic coverage goal; not forced)",
    )

    report = {
        "all_expected_matrix_checks_passed": all(
            check["passed"]
            for check in checks
            if check["check"] != "native triggered relocation was observed"
        ),
        "native_trigger_coverage_observed": triggered_count > 0,
        "note": (
            "The native-trigger check is reported as coverage information, not as "
            "a reason to alter the unconstrained physical simulation."
        ),
        "checks": checks,
    }
    with (output_dir / "validation_report.json").open("w") as stream:
        json.dump(report, stream, indent=2)
    with (output_dir / "validation_checks.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=["check", "passed", "observed", "expected"]
        )
        writer.writeheader()
        writer.writerows(checks)
    return checks


def group_legend() -> list[Line2D]:
    return [
        Line2D([0], [0], color=GROUP_COLORS[name], lw=3, label=name)
        for name, _ in GROUPS
    ]


def plot_overview(output_dir: Path, arrays: dict[str, np.ndarray]) -> None:
    years = arrays["year"]
    figure, axes = plt.subplots(3, 2, figsize=(17, 15), constrained_layout=True)
    for domain in range(DOMAIN_COUNT):
        group = GROUP_FOR_DOMAIN[domain]
        color = GROUP_COLORS[group]
        label = f"D{domain} ({ROAD_SETBACKS_M[domain]:g} m)"
        axes[0, 0].plot(
            years,
            arrays["shoreline_m"][domain] - arrays["shoreline_m"][domain, 0],
            color=color,
            alpha=0.72,
            lw=1.3,
            label=label,
        )
        axes[0, 1].plot(
            years, arrays["beach_width_m"][domain], color=color, alpha=0.72, lw=1.3
        )
        axes[1, 0].plot(
            years, arrays["road_setback_m"][domain], color=color, alpha=0.72, lw=1.3
        )
        axes[1, 1].plot(
            years,
            arrays["dune_crest_mean_m_MHW"][domain],
            color=color,
            alpha=0.72,
            lw=1.3,
        )
        axes[2, 0].plot(
            years,
            np.cumsum(np.nan_to_num(arrays["overwash_flux_m3_per_m"][domain])),
            color=color,
            alpha=0.72,
            lw=1.3,
        )

    axes[0, 0].axhline(0, color="black", lw=0.8)
    axes[0, 0].set(title="Shoreline position relative to year 0", ylabel="Change (m)")
    axes[0, 1].set(title="Roadway-managed beach width", ylabel="Beach width (m)")
    axes[1, 0].axhline(
        ROAD_SETBACK_TRIGGER_M,
        color="#cb181d",
        ls="--",
        lw=2,
        label="20 m dune-rebuild boundary",
    )
    axes[1, 0].set(title="Road setback", ylabel="Road-to-dune distance (m)")
    axes[1, 1].set(title="Mean dune crest elevation", ylabel="Elevation (m MHW)")
    axes[2, 0].set(
        title="Cumulative Barrier3D overwash flux", ylabel="Cumulative Qow (m³/m)"
    )

    event_axis = axes[2, 1]
    event_specs = (
        ("nourishment", "o", "#2171b5", "Nourishment"),
        ("dunes_rebuilt", "+", "#e31a1c", "Dunes rebuilt"),
        ("triggered_relocation", "*", "#ff7f00", "Triggered relocation"),
        ("forced_relocation", "s", "#33a02c", "Historical/forced relocation"),
        ("incomplete_relocation", "x", "#6a3d9a", "Incomplete relocation"),
    )
    for key, marker, color, label in event_specs:
        domain_indices, year_indices = np.where(arrays[key])
        event_axis.scatter(
            year_indices,
            domain_indices,
            marker=marker,
            color=color,
            s=45,
            linewidths=1.5,
            label=label,
        )
    event_axis.set(
        title="Calculated and requested management events",
        ylabel="Synthetic domain",
        yticks=np.arange(DOMAIN_COUNT),
        ylim=(-0.75, DOMAIN_COUNT - 0.25),
    )
    event_axis.legend(loc="upper left", bbox_to_anchor=(1.01, 1), fontsize=9)

    for axis in axes.flat:
        axis.set_xlabel("Model year")
        axis.grid(alpha=0.2)
    axes[0, 0].legend(
        handles=group_legend(), loc="upper left", fontsize=9, title="Domain group"
    )
    axes[1, 0].legend(
        handles=group_legend()
        + [Line2D([0], [0], color="#cb181d", ls="--", lw=2, label="20 m boundary")],
        loc="upper left",
        fontsize=8,
    )
    figure.suptitle(
        "20-domain full RoadwayManager synthetic integration test\n"
        "Active Barrier3D storms and RSLR; no prescribed physical trajectory",
        fontsize=17,
        fontweight="bold",
    )
    figure.savefig(output_dir / "full_functionality_overview.png", dpi=180)
    plt.close(figure)


def plot_heatmaps(output_dir: Path, arrays: dict[str, np.ndarray]) -> None:
    panels = (
        ("road_setback_m", "Road setback (m)", "viridis"),
        ("beach_width_m", "Beach width (m)", "YlGnBu"),
        ("dune_crest_mean_m_MHW", "Mean dune crest (m MHW)", "terrain"),
        ("mean_interior_height_m_MHW", "Mean interior elevation (m MHW)", "cividis"),
        ("shoreline_m", "Shoreline position (m)", "coolwarm"),
        ("overwash_flux_m3_per_m", "Annual overwash flux Qow (m³/m)", "magma"),
    )
    figure, axes = plt.subplots(3, 2, figsize=(17, 14), constrained_layout=True)
    for axis, (key, title, cmap) in zip(axes.flat, panels):
        values = arrays[key].copy()
        if key == "shoreline_m":
            values -= values[:, [0]]
            title = "Shoreline change from year 0 (m)"
        image = axis.imshow(values, aspect="auto", interpolation="nearest", cmap=cmap)
        axis.set(title=title, ylabel="Synthetic domain", xlabel="Model year")
        axis.set_yticks(np.arange(DOMAIN_COUNT))
        figure.colorbar(image, ax=axis, shrink=0.86)
        for boundary in (3.5, 7.5, 11.5, 15.5):
            axis.axhline(boundary, color="white", lw=1.2, alpha=0.9)
    figure.suptitle(
        "All 20 synthetic domains — annual calculated state", fontsize=17, fontweight="bold"
    )
    figure.savefig(output_dir / "all_domains_state_heatmaps.png", dpi=180)
    plt.close(figure)


def event_text(arrays: dict[str, np.ndarray], domain: int, year: int) -> str:
    labels = []
    if arrays["nourishment"][domain, year]:
        labels.append("N")
    if arrays["dunes_rebuilt"][domain, year]:
        labels.append("D")
    if arrays["triggered_relocation"][domain, year]:
        labels.append("T")
    if arrays["forced_relocation"][domain, year]:
        labels.append("H")
    if arrays["incomplete_relocation"][domain, year]:
        labels.append("X")
    return " ".join(labels)


def plot_domain_tiles(
    output_dir: Path,
    arrays: dict[str, np.ndarray],
    year: int,
    filename: str,
) -> None:
    figure, axes = plt.subplots(4, 5, figsize=(18, 11))
    figure.subplots_adjust(
        left=0.03, right=0.99, bottom=0.105, top=0.90, wspace=0.06, hspace=0.28
    )
    all_shore = arrays["shoreline_m"] - arrays["dune_grid_m"][:, [0]]
    all_dune = arrays["dune_grid_m"] - arrays["dune_grid_m"][:, [0]]
    road = all_dune + arrays["road_setback_m"]
    finite_road = road[np.isfinite(road)]
    x_min = float(np.nanmin(all_shore) - 15)
    x_max = float(max(np.nanmax(all_dune) + 60, np.nanmax(finite_road) + 20))

    for domain, axis in enumerate(axes.flat):
        shoreline = all_shore[domain, year]
        dune = all_dune[domain, year]
        road_position = road[domain, year]
        axis.axvspan(x_min, shoreline, color="#9ecae1", alpha=0.95)
        if shoreline < dune:
            axis.axvspan(shoreline, dune, color="#fdd49e", alpha=0.95)
        axis.axvspan(dune, x_max, color="#c7e9c0", alpha=0.9)
        axis.axvline(shoreline, color="#08519c", lw=2)
        axis.axvline(dune, color="#8c510a", lw=3)
        if np.isfinite(road_position):
            axis.axvspan(road_position, road_position + 3, color="#525252", alpha=0.95)
        rebuild_status = (
            "natural (>20 m)"
            if arrays["rebuild_disabled_by_setback"][domain, year]
            else "eligible (≤20 m)"
        )
        axis.text(
            0.02,
            0.95,
            f"D{domain}: {GROUP_FOR_DOMAIN[domain]}\n"
            f"setback={arrays['road_setback_m'][domain, year]:.1f} m | {rebuild_status}\n"
            f"beach={arrays['beach_width_m'][domain, year]:.1f} m | "
            f"crest={arrays['dune_crest_mean_m_MHW'][domain, year]:.2f} m MHW\n"
            f"storms={arrays['storm_count'][domain, year]} | "
            f"Qow={arrays['overwash_flux_m3_per_m'][domain, year]:.1f} | "
            f"events={event_text(arrays, domain, year) or '—'}",
            transform=axis.transAxes,
            va="top",
            fontsize=7.5,
            bbox={"facecolor": "white", "alpha": 0.88, "edgecolor": GROUP_COLORS[GROUP_FOR_DOMAIN[domain]]},
        )
        axis.set_xlim(x_min, x_max)
        axis.set_ylim(0, 1)
        axis.set_yticks([])
        axis.set_xlabel("Relative cross-shore position (m)", fontsize=7)
        axis.tick_params(axis="x", labelsize=7)

    figure.legend(
        handles=[
            Patch(color="#9ecae1", label="Ocean"),
            Patch(color="#fdd49e", label="Beach"),
            Patch(color="#c7e9c0", label="Barrier interior"),
            Line2D([0], [0], color="#8c510a", lw=3, label="Dune line"),
            Patch(color="#525252", label="Road"),
            Line2D([0], [0], marker="o", color="none", label="N=nourishment; D=dune rebuild; T=triggered relocation; H=historical/forced; X=incomplete"),
        ],
        loc="lower center",
        bbox_to_anchor=(0.5, 0.012),
        ncol=3,
        fontsize=9,
        frameon=True,
    )
    figure.suptitle(
        f"20 synthetic Barrier3D domains — model year {year}\n"
        "Calculated shoreline/dune positions and RoadwayManager state",
        fontsize=16,
        fontweight="bold",
    )
    figure.savefig(output_dir / filename, dpi=160)
    plt.close(figure)


def make_slow_gif(output_dir: Path, arrays: dict[str, np.ndarray]) -> None:
    final_year = int(arrays["year"][-1])
    event_years = set()
    for key in (
        "nourishment",
        "dunes_rebuilt",
        "triggered_relocation",
        "forced_relocation",
        "incomplete_relocation",
    ):
        event_years.update(np.where(arrays[key])[1].tolist())
    frames = sorted(set(range(0, final_year + 1, 5)) | event_years | {final_year})

    figure, axes = plt.subplots(4, 5, figsize=(18, 11), constrained_layout=True)
    all_shore = arrays["shoreline_m"] - arrays["dune_grid_m"][:, [0]]
    all_dune = arrays["dune_grid_m"] - arrays["dune_grid_m"][:, [0]]
    road = all_dune + arrays["road_setback_m"]
    finite_road = road[np.isfinite(road)]
    x_min = float(np.nanmin(all_shore) - 15)
    x_max = float(max(np.nanmax(all_dune) + 60, np.nanmax(finite_road) + 20))

    def update(frame_number: int):
        year = frames[frame_number]
        for domain, axis in enumerate(axes.flat):
            axis.clear()
            shoreline = all_shore[domain, year]
            dune = all_dune[domain, year]
            road_position = road[domain, year]
            axis.axvspan(x_min, shoreline, color="#9ecae1", alpha=0.95)
            if shoreline < dune:
                axis.axvspan(shoreline, dune, color="#fdd49e", alpha=0.95)
            axis.axvspan(dune, x_max, color="#c7e9c0", alpha=0.9)
            axis.axvline(shoreline, color="#08519c", lw=2)
            axis.axvline(dune, color="#8c510a", lw=3)
            if np.isfinite(road_position):
                axis.axvspan(
                    road_position, road_position + 3, color="#525252", alpha=0.95
                )
            events = event_text(arrays, domain, year)
            edge = "#e31a1c" if events else GROUP_COLORS[GROUP_FOR_DOMAIN[domain]]
            setback = arrays["road_setback_m"][domain, year]
            setback_label = f"{setback:.1f} m" if np.isfinite(setback) else "road stopped"
            axis.text(
                0.02,
                0.95,
                f"D{domain} | {GROUP_FOR_DOMAIN[domain]}\n"
                f"setback {setback_label}; beach {arrays['beach_width_m'][domain, year]:.1f} m\n"
                f"crest {arrays['dune_crest_mean_m_MHW'][domain, year]:.2f} m MHW; "
                f"storms {arrays['storm_count'][domain, year]}\n"
                f"events: {events or '—'}",
                transform=axis.transAxes,
                va="top",
                fontsize=7.5,
                bbox={"facecolor": "white", "alpha": 0.9, "edgecolor": edge, "linewidth": 1.5},
            )
            axis.set_xlim(x_min, x_max)
            axis.set_ylim(0, 1)
            axis.set_yticks([])
            axis.tick_params(axis="x", labelsize=7)
            axis.set_xlabel("Relative cross-shore position (m)", fontsize=7)
        figure.suptitle(
            f"20-domain full RoadwayManager test — model year {year}\n"
            "N=nourishment | D=dune rebuild | T=triggered relocation | "
            "H=historical/forced relocation | X=incomplete",
            fontsize=15,
            fontweight="bold",
        )
        return []

    animation = FuncAnimation(figure, update, frames=len(frames), blit=False)
    animation.save(
        output_dir / "all_domains_evolution_slow.gif",
        writer=PillowWriter(fps=1),
        dpi=90,
    )
    plt.close(figure)


def print_summary(arrays: dict[str, np.ndarray], cascade: Cascade, stop_reason: str) -> None:
    print("\nRun summary", flush=True)
    print(f"  final recorded year: {arrays['year'][-1]}", flush=True)
    print(f"  stop reason: {stop_reason}", flush=True)
    print(
        f"  active roads at end: {sum(not bool(v) for v in cascade.road_break)}/{DOMAIN_COUNT}",
        flush=True,
    )
    for group, domains_range in GROUPS:
        domains = list(domains_range)
        nourishments = int(np.sum(arrays["nourishment"][domains]))
        triggered = int(np.sum(arrays["triggered_relocation"][domains]))
        forced = int(np.sum(arrays["forced_relocation"][domains]))
        incomplete = int(np.sum(arrays["incomplete_relocation"][domains]))
        rebuilds = int(np.sum(arrays["dunes_rebuilt"][domains]))
        print(
            f"  {group:24s}: nourishment={nourishments}, "
            f"triggered relocation={triggered}, forced relocation={forced}, "
            f"incomplete={incomplete}, dune rebuild={rebuilds}",
            flush=True,
        )


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    input_dir = prepare_fresh_run(output_dir)
    validation = validate_inputs(input_dir)
    write_manifest(output_dir, input_dir, validation, status="running")

    cascade = None
    try:
        cascade, completed_year, stop_reason = run_model(input_dir)
        arrays = extract_results(cascade, completed_year)
        np.savez_compressed(output_dir / "cascade_run.npz", cascade=cascade)
        np.savez_compressed(output_dir / "annual_metrics.npz", **arrays)
        write_summary_csv(output_dir, arrays, cascade)
        write_validation_report(output_dir, arrays, cascade)
        plot_overview(output_dir, arrays)
        plot_heatmaps(output_dir, arrays)
        plot_domain_tiles(
            output_dir,
            arrays,
            int(arrays["year"][-1]),
            "all_domains_final_state.png",
        )
        if not args.skip_gif:
            print("Rendering slow GIF...", flush=True)
            make_slow_gif(output_dir, arrays)
        write_manifest(
            output_dir,
            input_dir,
            validation,
            status="complete",
            completed_year=int(arrays["year"][-1]),
            stop_reason=stop_reason,
        )
        print_summary(arrays, cascade, stop_reason)
        print(f"\nOutputs: {output_dir}", flush=True)
    except Exception as error:
        completed = None
        if cascade is not None:
            completed = min(len(barrier.x_s_TS) for barrier in cascade.barrier3d) - 1
        write_manifest(
            output_dir,
            input_dir,
            validation,
            status="failed",
            completed_year=completed,
            stop_reason=f"{type(error).__name__}: {error}",
        )
        raise


if __name__ == "__main__":
    main()
