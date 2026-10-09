#!/usr/bin/env python3
"""Run one 20-segment, alongshore-connected synthetic barrier island.

The simulation uses one CASCADE object with BRIE alongshore transport enabled.
Each Barrier3D segment represents an adjacent 500 m reach, giving a continuous
10 km synthetic island. Historical nourishment and relocation requests are
test inputs; the model calculates their feasibility and all physical responses.

This script never overwrites an existing output directory and does not modify
model source files.
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
from matplotlib.patches import Patch, Rectangle

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


SOURCE_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = SOURCE_ROOT.parent
sys.path.insert(0, str(SOURCE_ROOT))

from cascade import Cascade  # noqa: E402


DOMAIN_COUNT = 20
DOMAIN_LENGTH_M = 500
ISLAND_LENGTH_KM = DOMAIN_COUNT * DOMAIN_LENGTH_M / 1000
TIME_STEP_COUNT = 200

ROAD_ELEVATION_M_MHW = 1.7
ROAD_WIDTH_M = 30.0
ROAD_SETBACK_M = 30.0
ROAD_SETBACK_TRIGGER_M = 20.0
DUNE_DESIGN_ELEVATION_M_MHW = 3.7
DUNE_MINIMUM_ELEVATION_M_MHW = 2.2
NOURISHMENT_VOLUME_M3_PER_M = 300.0

# A transparent synthetic management history. Domain ranges are inclusive.
SYNTHETIC_NOURISHMENT = {
    20: (0, 4),
    60: (5, 9),
    100: (10, 14),
    140: (15, 19),
    180: (0, 19),
}
SYNTHETIC_RELOCATION = {
    40: (3, 6),
    80: (8, 11),
    120: (13, 16),
    160: (17, 19),
}
HISTORICAL_RELOCATION_TARGET_M = 30.0

DATA_ROOT = SOURCE_ROOT / "data"
INPUT_SOURCES = {
    "barrier3d-default-elevation.npy": DATA_ROOT / "barrier3d-default-elevation.npy",
    "barrier3d-default-dunes.npy": DATA_ROOT / "barrier3d-default-dunes.npy",
    "barrier3d-default-parameters.yaml": DATA_ROOT / "barrier3d-default-parameters.yaml",
    "cascade-default-storms.npy": DATA_ROOT / "cascade-default-storms.npy",
    # The default YAML requires this path during validation. CASCADE's
    # initialize_equal subsequently uses default rmin/rmax initialization.
    "barrier3d-default-growthparam.npy": SOURCE_ROOT
    / "tests"
    / "test_human_dynamics"
    / "growthparam_1000dam.npy",
}
INPUT_FILES = tuple(INPUT_SOURCES)
OUTPUT_DIR = (
    WORKSPACE_ROOT
    / "CASCADE"
    / "output"
    / "roadway_functionalities_full_synthetic"
    / "connected_20_domain_island_200_years"
)
GIF_FPS = 1.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--road-width",
        type=float,
        default=ROAD_WIDTH_M,
        help="Modeled road width in meters (default: 30).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=OUTPUT_DIR,
        help="Fresh output directory; existing directories are rejected.",
    )
    parser.add_argument(
        "--gif-fps",
        type=float,
        default=GIF_FPS,
        help="GIF frames per second; 0.5 gives two seconds per frame.",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def expanded_domains(bounds: tuple[int, int]) -> list[int]:
    return list(range(bounds[0], bounds[1] + 1))


def prepare_connected_inputs() -> tuple[Path, list[str], list[str], dict]:
    if OUTPUT_DIR.exists():
        raise FileExistsError(f"Refusing to overwrite existing run: {OUTPUT_DIR}")
    input_dir = OUTPUT_DIR / "inputs"
    input_dir.mkdir(parents=True)
    for filename, source in INPUT_SOURCES.items():
        if not source.is_file():
            raise FileNotFoundError(source)
        shutil.copy2(source, input_dir / filename)

    elevation = np.load(input_dir / "barrier3d-default-elevation.npy")
    dunes = np.load(input_dir / "barrier3d-default-dunes.npy")
    growth = np.load(input_dir / "barrier3d-default-growthparam.npy")
    if elevation.ndim != 2 or elevation.shape[1] < 50:
        raise ValueError(f"Unexpected default elevation shape: {elevation.shape}")
    if dunes.size < DOMAIN_COUNT * 50:
        raise ValueError(
            f"Need {DOMAIN_COUNT * 50} dune values; found {dunes.size}."
        )

    # Build one 1,000-cell alongshore synthetic initial field. The default
    # elevation is tiled only to provide sufficient alongshore coverage, then
    # both elevation and dune fields are split into adjacent, non-overlapping
    # 50-cell reaches for the 20 Barrier3D models.
    repeats = int(np.ceil((DOMAIN_COUNT * 50) / elevation.shape[1]))
    connected_elevation = np.tile(elevation, (1, repeats))[:, : DOMAIN_COUNT * 50]
    elevation_files = []
    dune_files = []
    for domain in range(DOMAIN_COUNT):
        start = domain * 50
        stop = start + 50
        elevation_name = f"connected_elevation_domain_{domain:02d}.npy"
        dune_name = f"connected_dunes_domain_{domain:02d}.npy"
        np.save(input_dir / elevation_name, connected_elevation[:, start:stop])
        # The current default YAML has DuneParamMultipleRows enabled and a
        # two-cell dune width. Expand the standard one-row values across both
        # rows without altering their elevation.
        dune_segment = np.repeat(dunes[start:stop, None], 2, axis=1)
        np.save(input_dir / dune_name, dune_segment)
        elevation_files.append(elevation_name)
        dune_files.append(dune_name)

    with (input_dir / "barrier3d-default-parameters.yaml").open() as stream:
        parameters = yaml.safe_load(stream)
    validation = {
        "base_elevation_shape": list(elevation.shape),
        "connected_elevation_shape": list(connected_elevation.shape),
        "base_dune_shape": list(dunes.shape),
        "growth_parameter_shape": list(growth.shape),
        "growth_parameter_loading_after_initialize_equal": False,
        "segment_elevation_shape": [45, 50],
        "segment_dune_shape": [50, 2],
        "storm_start": int(parameters["StormStart"]),
        "parameter_TMAX_before_CASCADE_initialization": int(parameters["TMAX"]),
    }
    if validation["storm_start"] >= TIME_STEP_COUNT:
        raise ValueError("The default parameter file would disable storms.")
    return input_dir, elevation_files, dune_files, validation


def cascade_arguments(elevation_files: list[str], dune_files: list[str]) -> dict:
    return {
        "elevation_file": elevation_files,
        "dune_file": dune_files,
        "parameter_file": "barrier3d-default-parameters.yaml",
        "storm_file": "cascade-default-storms.npy",
        "num_cores": 1,
        "roadway_management_module": [True] * DOMAIN_COUNT,
        "alongshore_transport_module": True,
        "beach_nourishment_module": [False] * DOMAIN_COUNT,
        "community_economics_module": False,
        "outwash_module": [False] * DOMAIN_COUNT,
        "alongshore_section_count": DOMAIN_COUNT,
        "time_step_count": TIME_STEP_COUNT,
        "road_ele": [ROAD_ELEVATION_M_MHW] * DOMAIN_COUNT,
        "road_width": [ROAD_WIDTH_M] * DOMAIN_COUNT,
        "road_setback": [ROAD_SETBACK_M] * DOMAIN_COUNT,
        "road_setback_trigger": [ROAD_SETBACK_TRIGGER_M] * DOMAIN_COUNT,
        "dune_design_elevation": [DUNE_DESIGN_ELEVATION_M_MHW] * DOMAIN_COUNT,
        "dune_minimum_elevation": [DUNE_MINIMUM_ELEVATION_M_MHW] * DOMAIN_COUNT,
        "nourishment_interval": [None] * DOMAIN_COUNT,
        "nourishment_volume": [NOURISHMENT_VOLUME_M3_PER_M] * DOMAIN_COUNT,
        "group_roadway_abandonment": None,
        "trigger_dune_knockdown": False,
    }


def write_manifest(
    input_dir: Path,
    elevation_files: list[str],
    dune_files: list[str],
    input_validation: dict,
    status: str,
    completed_year: int | None = None,
    stop_reason: str | None = None,
) -> None:
    manifest = {
        "status": status,
        "description": "one 10-km synthetic island with 20 connected 500-m segments",
        "source_root": str(SOURCE_ROOT),
        "roadway_manager_source": str(SOURCE_ROOT / "cascade" / "roadway_manager.py"),
        "domain_count": DOMAIN_COUNT,
        "domain_length_m": DOMAIN_LENGTH_M,
        "island_length_km": ISLAND_LENGTH_KM,
        "requested_states": TIME_STEP_COUNT,
        "gif_frames_per_second": GIF_FPS,
        "completed_year": completed_year,
        "stop_reason": stop_reason,
        "connectivity": {
            "one_Cascade_object": True,
            "BRIE_alongshore_transport": True,
            "Barrier3D_segments_are_adjacent": True,
        },
        "management": {
            "RoadwayManager": "enabled on all domains",
            "BeachDuneManager": "disabled on all domains",
            "Outwasher": "disabled on all domains",
            "automatic_nourishment_interval": None,
            "historical_nourishment_volume_m3_per_m": NOURISHMENT_VOLUME_M3_PER_M,
            "historical_relocation_target_setback_m": HISTORICAL_RELOCATION_TARGET_M,
            "native_triggered_relocation": "enabled by RoadwayManager",
            "road_dune_rebuild_rule": "eligible at setback <= 20 m; skipped at > 20 m",
        },
        "synthetic_historical_nourishment": {
            str(year): expanded_domains(bounds)
            for year, bounds in SYNTHETIC_NOURISHMENT.items()
        },
        "synthetic_historical_relocation": {
            str(year): expanded_domains(bounds)
            for year, bounds in SYNTHETIC_RELOCATION.items()
        },
        "initial_conditions": {
            "road_elevation_m_MHW": ROAD_ELEVATION_M_MHW,
            "road_width_m": ROAD_WIDTH_M,
            "road_setback_m": ROAD_SETBACK_M,
            "road_setback_trigger_m": ROAD_SETBACK_TRIGGER_M,
            "dune_design_elevation_m_MHW": DUNE_DESIGN_ELEVATION_M_MHW,
            "dune_minimum_elevation_m_MHW": DUNE_MINIMUM_ELEVATION_M_MHW,
            "RSLR_m_per_year": 0.004,
            "RSLR_constant": True,
            "background_erosion_m_per_year": 0.0,
            "wave_height_m": 1,
            "wave_period_s": 7,
            "MHW_m_NAVD88": 0.46,
            "berm_elevation_m_NAVD88": 1.9,
            "beach_slope": 0.04,
        },
        "input_validation": input_validation,
        "segment_elevation_files": elevation_files,
        "segment_dune_files": dune_files,
        "copied_input_sha256": {
            filename: sha256(input_dir / filename) for filename in INPUT_FILES
        },
        "no_prescribed_physical_state": True,
    }
    with (OUTPUT_DIR / "manifest.json").open("w") as stream:
        json.dump(manifest, stream, indent=2)


def run_model(input_dir: Path, elevation_files: list[str], dune_files: list[str]):
    cascade = Cascade(
        str(input_dir),
        name="connected_20_domain_synthetic_island",
        **cascade_arguments(elevation_files, dune_files),
    )
    event_log = []
    completed_year = 0
    stop_reason = "completed requested duration"

    for year in range(1, TIME_STEP_COUNT):
        nourish_domains = expanded_domains(SYNTHETIC_NOURISHMENT[year]) if year in SYNTHETIC_NOURISHMENT else []
        relocation_domains = expanded_domains(SYNTHETIC_RELOCATION[year]) if year in SYNTHETIC_RELOCATION else []

        cascade.nourish_now = [
            int(domain in nourish_domains) for domain in range(DOMAIN_COUNT)
        ]
        for domain in relocation_domains:
            if not cascade.road_break[domain]:
                cascade.roadways[domain].request_historical_relocation(
                    HISTORICAL_RELOCATION_TARGET_M
                )
                event_log.append(
                    {
                        "year": year,
                        "domain": domain,
                        "request": "historical_relocation",
                        "requested_value": HISTORICAL_RELOCATION_TARGET_M,
                        "unit": "m setback",
                        "queued": True,
                    }
                )
            else:
                event_log.append(
                    {
                        "year": year,
                        "domain": domain,
                        "request": "historical_relocation",
                        "requested_value": HISTORICAL_RELOCATION_TARGET_M,
                        "unit": "m setback",
                        "queued": False,
                    }
                )
        for domain in nourish_domains:
            event_log.append(
                {
                    "year": year,
                    "domain": domain,
                    "request": "historical_nourishment",
                    "requested_value": NOURISHMENT_VOLUME_M3_PER_M,
                    "unit": "m3/m",
                    "queued": not bool(cascade.road_break[domain]),
                }
            )

        cascade.update()
        completed_year = min(len(b.x_s_TS) for b in cascade.barrier3d) - 1
        if cascade.b3d_break:
            drowned = [i for i, b in enumerate(cascade.barrier3d) if b.drown_break]
            stop_reason = f"Barrier3D drowning; domains={drowned}"
            print(f"Stopped at year {completed_year}: {stop_reason}", flush=True)
            break
        if year == 1 or year % 10 == 0 or year == TIME_STEP_COUNT - 1:
            active = sum(not bool(value) for value in cascade.road_break)
            print(
                f"Completed year {year:3d}/{TIME_STEP_COUNT - 1}; active roads={active}/20",
                flush=True,
            )

    return cascade, completed_year, stop_reason, event_log


def pad(values, count: int, fill=np.nan, dtype=float) -> np.ndarray:
    result = np.full(count, fill, dtype=dtype)
    array = np.asarray(values)
    copied = min(count, array.size)
    result[:copied] = array[:copied]
    return result


def extract_metrics(cascade: Cascade, completed_year: int) -> dict[str, np.ndarray]:
    count = completed_year + 1
    shape = (DOMAIN_COUNT, count)
    metrics = {
        "shoreline_m": np.full(shape, np.nan),
        "shoreface_toe_m": np.full(shape, np.nan),
        "beach_width_m": np.full(shape, np.nan),
        "road_setback_m": np.full(shape, np.nan),
        "road_elevation_m_MHW": np.full(shape, np.nan),
        "dune_crest_mean_m_MHW": np.full(shape, np.nan),
        "interior_height_m_MHW": np.full(shape, np.nan),
        "storm_count": np.zeros(shape, dtype=int),
        "overwash_flux_m3_per_m": np.zeros(shape),
        "nourishment": np.zeros(shape, dtype=bool),
        "nourishment_volume_m3_per_m": np.zeros(shape),
        "triggered_relocation": np.zeros(shape, dtype=bool),
        "historical_request": np.zeros(shape, dtype=bool),
        "forced_relocation": np.zeros(shape, dtype=bool),
        "incomplete_relocation": np.zeros(shape, dtype=bool),
        "dunes_rebuilt": np.zeros(shape, dtype=bool),
        "rebuild_disabled_by_setback": np.zeros(shape, dtype=bool),
        "dune_migration_on": np.ones(shape, dtype=bool),
    }
    for domain, (barrier, manager) in enumerate(zip(cascade.barrier3d, cascade.roadways)):
        metrics["shoreline_m"][domain] = pad(np.asarray(barrier.x_s_TS) * 10, count)
        metrics["shoreface_toe_m"][domain] = pad(np.asarray(barrier.x_t_TS) * 10, count)
        metrics["beach_width_m"][domain] = pad(manager.beach_width, count)
        road_width = pad(manager._road_width_TS, count)
        setback = pad(manager._road_setback_TS, count)
        road_elevation = pad(manager._road_ele_TS, count)
        setback[road_width <= 0] = np.nan
        road_elevation[road_width <= 0] = np.nan
        metrics["road_setback_m"][domain] = setback
        metrics["road_elevation_m_MHW"][domain] = road_elevation
        dunes = np.asarray(barrier.DuneDomain[:count])
        crest = np.mean(np.max(dunes, axis=2), axis=1) * 10 + barrier.BermEl * 10
        metrics["dune_crest_mean_m_MHW"][domain, : crest.size] = crest
        metrics["interior_height_m_MHW"][domain] = pad(
            np.asarray(barrier.h_b_TS) * 10, count
        )
        metrics["storm_count"][domain] = pad(
            barrier._StormCount, count, fill=0, dtype=int
        )
        metrics["overwash_flux_m3_per_m"][domain] = pad(
            barrier.QowTS, count, fill=0
        )
        for key, values in (
            ("nourishment", manager.nourishment_TS),
            ("triggered_relocation", manager.triggered_relocation_TS),
            ("historical_request", manager.historical_relocation_requested_TS),
            ("forced_relocation", manager.forced_relocation_TS),
            ("incomplete_relocation", manager.relocation_incomplete_TS),
            ("dunes_rebuilt", manager._dunes_rebuilt_TS),
            ("rebuild_disabled_by_setback", manager.road_dune_rebuild_disabled_TS),
        ):
            metrics[key][domain] = pad(values, count, fill=False, dtype=bool)
        metrics["nourishment_volume_m3_per_m"][domain] = pad(
            manager.nourishment_volume_TS, count, fill=0
        )
        migration = pad(manager.dune_migration_on, count)
        metrics["dune_migration_on"][domain] = np.nan_to_num(migration, nan=1).astype(bool)
    metrics["year"] = np.arange(count)
    metrics["domain"] = np.arange(DOMAIN_COUNT)
    return metrics


def write_event_log(event_log: list[dict]) -> None:
    with (OUTPUT_DIR / "synthetic_historical_requests.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=["year", "domain", "request", "requested_value", "unit", "queued"],
        )
        writer.writeheader()
        writer.writerows(event_log)


def years(values: np.ndarray) -> str:
    return ";".join(str(value) for value in np.flatnonzero(values))


def write_domain_summary(metrics: dict[str, np.ndarray], cascade: Cascade) -> None:
    fields = [
        "domain", "alongshore_start_km", "alongshore_stop_km", "road_active_at_end",
        "nourishment_years", "historical_request_years", "forced_relocation_years",
        "triggered_relocation_years", "incomplete_relocation_years", "dune_rebuild_years",
        "final_road_setback_m", "final_beach_width_m", "final_dune_crest_m_MHW",
        "shoreline_change_m", "total_storms", "cumulative_overwash_flux_m3_per_m",
    ]
    with (OUTPUT_DIR / "domain_summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for domain in range(DOMAIN_COUNT):
            shoreline = metrics["shoreline_m"][domain]
            writer.writerow(
                {
                    "domain": domain,
                    "alongshore_start_km": domain * 0.5,
                    "alongshore_stop_km": (domain + 1) * 0.5,
                    "road_active_at_end": not bool(cascade.road_break[domain]),
                    "nourishment_years": years(metrics["nourishment"][domain]),
                    "historical_request_years": years(metrics["historical_request"][domain]),
                    "forced_relocation_years": years(metrics["forced_relocation"][domain]),
                    "triggered_relocation_years": years(metrics["triggered_relocation"][domain]),
                    "incomplete_relocation_years": years(metrics["incomplete_relocation"][domain]),
                    "dune_rebuild_years": years(metrics["dunes_rebuilt"][domain]),
                    "final_road_setback_m": metrics["road_setback_m"][domain, -1],
                    "final_beach_width_m": metrics["beach_width_m"][domain, -1],
                    "final_dune_crest_m_MHW": metrics["dune_crest_mean_m_MHW"][domain, -1],
                    "shoreline_change_m": shoreline[-1] - shoreline[0],
                    "total_storms": int(np.sum(metrics["storm_count"][domain])),
                    "cumulative_overwash_flux_m3_per_m": float(
                        np.sum(metrics["overwash_flux_m3_per_m"][domain])
                    ),
                }
            )


def validate_run(metrics: dict[str, np.ndarray], cascade: Cascade) -> dict:
    checks = []
    def add(name: str, passed: bool, observed, expected) -> None:
        checks.append(
            {"check": name, "passed": bool(passed), "observed": str(observed), "expected": str(expected)}
        )

    add("200 annual states saved", metrics["year"].size == 200, metrics["year"].size, 200)
    add("one Cascade contains 20 Barrier3D segments", len(cascade.barrier3d) == 20, len(cascade.barrier3d), 20)
    add(
        "BRIE alongshore transport enabled",
        bool(cascade._alongshore_transport_module),
        cascade._alongshore_transport_module,
        True,
    )
    add("storms active", int(np.sum(metrics["storm_count"])) > 0, int(np.sum(metrics["storm_count"])), ">0")

    requested_nourishment = sum(
        len(expanded_domains(bounds)) for bounds in SYNTHETIC_NOURISHMENT.values()
    )
    observed_nourishment = int(np.sum(metrics["nourishment"]))
    add("historical nourishment requests applied", observed_nourishment == requested_nourishment, observed_nourishment, requested_nourishment)
    nourishment_response = []
    for domain, year_index in zip(*np.where(metrics["nourishment"])):
        nourishment_response.append(
            metrics["shoreline_m"][domain, year_index]
            < metrics["shoreline_m"][domain, year_index - 1]
        )
    add(
        "nourishment events produced net same-year progradation",
        all(nourishment_response),
        f"{sum(nourishment_response)}/{len(nourishment_response)}",
        f"{len(nourishment_response)}/{len(nourishment_response)}",
    )

    requested_relocation = sum(
        len(expanded_domains(bounds)) for bounds in SYNTHETIC_RELOCATION.values()
    )
    observed_requests = int(np.sum(metrics["historical_request"]))
    observed_forced = int(np.sum(metrics["forced_relocation"]))
    add("historical relocation requests recorded", observed_requests == requested_relocation, observed_requests, requested_relocation)
    add("historical relocation requests succeeded", observed_forced == requested_relocation, observed_forced, requested_relocation)
    add("no incomplete relocation", not np.any(metrics["incomplete_relocation"]), int(np.sum(metrics["incomplete_relocation"])), 0)
    add("all roads active at final year", not any(cascade.road_break), sum(not bool(v) for v in cascade.road_break), 20)

    report = {
        "passed": all(check["passed"] for check in checks),
        "native_triggered_relocations_observed": int(np.sum(metrics["triggered_relocation"])),
        "dune_rebuild_events_observed": int(np.sum(metrics["dunes_rebuilt"])),
        "checks": checks,
    }
    with (OUTPUT_DIR / "validation_report.json").open("w") as stream:
        json.dump(report, stream, indent=2)
    return report


def plot_connected_metrics(metrics: dict[str, np.ndarray]) -> None:
    figure, axes = plt.subplots(3, 2, figsize=(17, 14), constrained_layout=True)
    panels = (
        ("shoreline_m", "Shoreline change from year 0 (m)", "coolwarm"),
        ("road_setback_m", "Road setback (m)", "viridis"),
        ("beach_width_m", "Beach width (m)", "YlGnBu"),
        ("dune_crest_mean_m_MHW", "Mean dune crest (m MHW)", "terrain"),
        ("interior_height_m_MHW", "Mean interior elevation (m MHW)", "cividis"),
        ("overwash_flux_m3_per_m", "Annual overwash flux Qow (m³/m)", "magma"),
    )
    for axis, (key, title, cmap) in zip(axes.flat, panels):
        values = metrics[key].copy()
        if key == "shoreline_m":
            values -= values[:, [0]]
        image = axis.imshow(
            values,
            origin="lower",
            aspect="auto",
            interpolation="nearest",
            extent=[0, metrics["year"][-1], 0, ISLAND_LENGTH_KM],
            cmap=cmap,
        )
        axis.set(title=title, xlabel="Model year", ylabel="Alongshore distance (km)")
        figure.colorbar(image, ax=axis, shrink=0.86)
        for year, bounds in SYNTHETIC_NOURISHMENT.items():
            axis.plot(year, (np.mean(bounds) + 0.5) * 0.5, "o", color="#00ffff", ms=4)
        for year, bounds in SYNTHETIC_RELOCATION.items():
            axis.plot(year, (np.mean(bounds) + 0.5) * 0.5, "*", color="#ff00ff", ms=7)
    figure.suptitle(
        "One connected 10-km synthetic barrier island\n"
        "Cyan circles: historical nourishment reaches | Magenta stars: historical relocation reaches",
        fontsize=16,
        fontweight="bold",
    )
    figure.savefig(OUTPUT_DIR / "connected_island_metrics.png", dpi=180)
    plt.close(figure)


def plan_view_limits(cascade: Cascade, metrics: dict[str, np.ndarray]) -> tuple[int, int]:
    max_width = max(
        domain.shape[0]
        for barrier in cascade.barrier3d
        for domain in barrier.DomainTS[: metrics["year"].size]
        if domain is not None
    )
    shoreline_cells = metrics["shoreline_m"] / 10
    beach_cells = np.nan_to_num(metrics["beach_width_m"], nan=30.0) / 10
    y_min = int(np.floor(np.nanmin(shoreline_cells))) - 8
    y_max = int(np.ceil(np.nanmax(shoreline_cells + beach_cells))) + max_width + 12
    return y_min, y_max


def connected_plan_view(
    cascade: Cascade,
    metrics: dict[str, np.ndarray],
    year: int,
    y_limits: tuple[int, int],
) -> tuple[np.ndarray, list[tuple[int, float, float]]]:
    y_min, y_max = y_limits
    barrier_length_cells = cascade.barrier3d[0].BarrierLength
    grid = np.full((y_max - y_min, barrier_length_cells * DOMAIN_COUNT), -1.0)
    road_rectangles = []

    for domain, (barrier, manager) in enumerate(zip(cascade.barrier3d, cascade.roadways)):
        shoreline = float(metrics["shoreline_m"][domain, year] / 10)
        beach_width = float(metrics["beach_width_m"][domain, year] / 10)
        if not np.isfinite(beach_width):
            beach_width = cascade._initial_beach_width[domain] / 10
        shore_cell = int(np.floor(shoreline))
        dune_toe = int(np.floor(shoreline + beach_width))
        cellular_beach_width = max(0, dune_toe - shore_cell)
        beach = np.zeros((cellular_beach_width, barrier_length_cells))
        if cellular_beach_width:
            increment = (barrier.BermEl - barrier.SL) / (cellular_beach_width + 1)
            for row in range(cellular_beach_width):
                beach[row] = (barrier.SL + increment * (row + 1)) * 10

        dunes = (barrier.DuneDomain[year] + barrier.BermEl) * 10
        dunes = np.flipud(np.rot90(dunes))
        interior = np.asarray(barrier.DomainTS[year]) * 10
        segment = np.vstack((beach, dunes, interior))
        segment[segment < 0] = -1
        row_start = shore_cell - y_min
        row_stop = min(grid.shape[0], row_start + segment.shape[0])
        segment = segment[: row_stop - row_start]
        x_start = domain * barrier_length_cells
        grid[row_start:row_stop, x_start : x_start + barrier_length_cells] = segment

        setback = metrics["road_setback_m"][domain, year]
        if np.isfinite(setback):
            road_y = dune_toe - y_min + setback / 10
            road_rectangles.append((x_start, road_y, ROAD_WIDTH_M / 10))
    return grid, road_rectangles


def compact_domain_ranges(domains) -> str:
    values = sorted(int(value) for value in domains)
    if not values:
        return "none"
    ranges = []
    start = previous = values[0]
    for value in values[1:]:
        if value == previous + 1:
            previous = value
            continue
        ranges.append(f"D{start}" if start == previous else f"D{start}–{previous}")
        start = previous = value
    ranges.append(f"D{start}" if start == previous else f"D{start}–{previous}")
    return ", ".join(ranges)


EVENT_OVERLAYS = (
    ("nourishment", "#00c8ff", "Nourishment applied"),
    ("forced_relocation", "#ff00c8", "Historical relocation completed"),
    ("triggered_relocation", "#ff8c00", "Native relocation completed"),
    ("incomplete_relocation", "#e31a1c", "Incomplete relocation"),
    ("dunes_rebuilt", "#7fff00", "Dunes rebuilt"),
)


def contiguous_runs(domains) -> list[tuple[int, int]]:
    values = sorted(int(value) for value in domains)
    if not values:
        return []
    runs = []
    start = previous = values[0]
    for value in values[1:]:
        if value == previous + 1:
            previous = value
        else:
            runs.append((start, previous))
            start = previous = value
    runs.append((start, previous))
    return runs


def add_event_reach_overlays(axis, metrics: dict[str, np.ndarray], year: int) -> None:
    for key, color, _ in EVENT_OVERLAYS:
        domains = np.flatnonzero(metrics[key][:, year])
        for start, stop in contiguous_runs(domains):
            axis.axvspan(
                start * 50,
                (stop + 1) * 50,
                facecolor=color,
                edgecolor=color,
                alpha=0.20,
                linewidth=3,
                zorder=5,
            )


def event_description(metrics: dict[str, np.ndarray], year: int) -> str:
    descriptions = []
    event_types = (
        ("nourishment", "nourishment"),
        ("forced_relocation", "historical relocation"),
        ("triggered_relocation", "native relocation"),
        ("incomplete_relocation", "incomplete relocation"),
        ("dunes_rebuilt", "dune rebuild"),
    )
    for key, label in event_types:
        domains = np.flatnonzero(metrics[key][:, year])
        if domains.size:
            descriptions.append(f"{label}: {compact_domain_ranges(domains)}")
    return " | ".join(descriptions) if descriptions else "No management event this year"


def gif_annotation_columns(metrics: dict[str, np.ndarray], year: int) -> tuple[str, str, str]:
    scheduled_nourishment = (
        expanded_domains(SYNTHETIC_NOURISHMENT[year])
        if year in SYNTHETIC_NOURISHMENT
        else []
    )
    scheduled_relocation = (
        expanded_domains(SYNTHETIC_RELOCATION[year])
        if year in SYNTHETIC_RELOCATION
        else []
    )
    applied_nourishment = np.flatnonzero(metrics["nourishment"][:, year])
    forced = np.flatnonzero(metrics["forced_relocation"][:, year])
    triggered = np.flatnonzero(metrics["triggered_relocation"][:, year])
    incomplete = np.flatnonzero(metrics["incomplete_relocation"][:, year])
    rebuilt = np.flatnonzero(metrics["dunes_rebuilt"][:, year])
    eligible = np.flatnonzero(
        np.isfinite(metrics["road_setback_m"][:, year])
        & (metrics["road_setback_m"][:, year] <= ROAD_SETBACK_TRIGGER_M)
    )
    natural = np.flatnonzero(metrics["rebuild_disabled_by_setback"][:, year])
    migration_off = np.flatnonzero(~metrics["dune_migration_on"][:, year])
    active = np.flatnonzero(np.isfinite(metrics["road_setback_m"][:, year]))

    requested = (
        "SYNTHETIC HISTORICAL INPUTS\n"
        f"Nourishment request: {compact_domain_ranges(scheduled_nourishment)}"
        + (
            f" at {NOURISHMENT_VOLUME_M3_PER_M:g} m³/m"
            if scheduled_nourishment
            else ""
        )
        + "\n"
        f"Relocation request: {compact_domain_ranges(scheduled_relocation)}"
        + (
            f" to {HISTORICAL_RELOCATION_TARGET_M:g} m setback"
            if scheduled_relocation
            else ""
        )
    )
    outcomes = (
        "MODEL-CALCULATED MANAGEMENT OUTPUTS\n"
        f"Nourishment applied: {compact_domain_ranges(applied_nourishment)}\n"
        f"Historical relocation completed: {compact_domain_ranges(forced)}\n"
        f"Native relocation completed: {compact_domain_ranges(triggered)}\n"
        f"Incomplete relocation: {compact_domain_ranges(incomplete)}\n"
        f"Dunes rebuilt: {compact_domain_ranges(rebuilt)}"
    )
    finite_beach = metrics["beach_width_m"][:, year]
    finite_beach = finite_beach[np.isfinite(finite_beach)]
    state = (
        "CURRENT FORCING AND MANAGEMENT STATE\n"
        f"Active roads: {len(active)}/20 | road width: {ROAD_WIDTH_M:g} m\n"
        f"Setback ≤20 m, rebuild eligible: {compact_domain_ranges(eligible)}\n"
        f"Setback >20 m, natural dunes: {compact_domain_ranges(natural)}\n"
        f"Dune migration OFF: {compact_domain_ranges(migration_off)}\n"
        f"Beach width range: {np.min(finite_beach):.1f}–{np.max(finite_beach):.1f} m\n"
        f"RSLR since year 0: {year * 0.004:.3f} m"
    )
    return requested, outcomes, state


def draw_connected_plan_view(
    axis,
    cascade: Cascade,
    metrics: dict[str, np.ndarray],
    year: int,
    y_limits: tuple[int, int],
) -> None:
    grid, roads = connected_plan_view(cascade, metrics, year, y_limits)
    axis.imshow(grid, origin="lower", aspect="auto", cmap="terrain", vmin=-1.1, vmax=4.0)
    add_event_reach_overlays(axis, metrics, year)
    for x_start, road_y, road_width in roads:
        axis.add_patch(
            Rectangle(
                (x_start, road_y),
                50,
                road_width,
                facecolor="#525252",
                edgecolor="white",
                linewidth=0.25,
                alpha=0.75,
            )
        )
    for domain in range(1, DOMAIN_COUNT):
        axis.axvline(domain * 50, color="white", lw=0.35, alpha=0.65)
    x_ticks = np.arange(0, DOMAIN_COUNT * 50 + 1, 100)
    axis.set_xticks(x_ticks, x_ticks / 100)
    y_min, y_max = y_limits
    y_ticks = np.linspace(0, y_max - y_min, 6)
    axis.set_yticks(y_ticks, np.round((y_ticks + y_min) * 10).astype(int))
    axis.set_xlabel("Alongshore distance (km)")
    axis.set_ylabel("Cross-shore coordinate (m)")
    axis.set_title(
        f"Model year {year} — {event_description(metrics, year)}",
        fontsize=11,
        fontweight="bold",
    )


def plot_connected_states(cascade: Cascade, metrics: dict[str, np.ndarray]) -> None:
    y_limits = plan_view_limits(cascade, metrics)
    figure, axes = plt.subplots(2, 1, figsize=(18, 8), constrained_layout=True)
    draw_connected_plan_view(axes[0], cascade, metrics, 0, y_limits)
    draw_connected_plan_view(axes[1], cascade, metrics, int(metrics["year"][-1]), y_limits)
    figure.suptitle(
        "Actual connected Barrier3D elevation domains and modeled roadway",
        fontsize=16,
        fontweight="bold",
    )
    figure.savefig(OUTPUT_DIR / "connected_island_initial_final.png", dpi=180)
    plt.close(figure)


def make_connected_gif(cascade: Cascade, metrics: dict[str, np.ndarray]) -> None:
    y_limits = plan_view_limits(cascade, metrics)
    event_years = set(SYNTHETIC_NOURISHMENT) | set(SYNTHETIC_RELOCATION)
    for key in ("triggered_relocation", "incomplete_relocation", "dunes_rebuilt"):
        event_years.update(np.where(metrics[key])[1].tolist())
    final_year = int(metrics["year"][-1])
    frames = sorted(set(range(0, final_year + 1, 5)) | event_years | {final_year})
    figure, (axis, annotation_axis) = plt.subplots(
        2,
        1,
        figsize=(18, 7.5),
        gridspec_kw={"height_ratios": [5.0, 1.75]},
        constrained_layout=True,
    )

    def update(frame_index: int):
        axis.clear()
        annotation_axis.clear()
        year = frames[frame_index]
        draw_connected_plan_view(axis, cascade, metrics, year, y_limits)
        storms = int(np.sum(metrics["storm_count"][:, year]))
        qow = float(np.sum(metrics["overwash_flux_m3_per_m"][:, year]))
        active = int(np.sum(np.isfinite(metrics["road_setback_m"][:, year])))
        axis.text(
            0.01,
            0.02,
            f"20 connected segments | active roads {active}/20 | storms {storms} | total Qow {qow:.1f} m³/m",
            transform=axis.transAxes,
            fontsize=9,
            bbox={"facecolor": "white", "alpha": 0.88, "edgecolor": "#333333"},
        )
        requested, outcomes, state = gif_annotation_columns(metrics, year)
        annotation_axis.axis("off")
        annotation_axis.text(
            0.01,
            0.98,
            requested,
            transform=annotation_axis.transAxes,
            va="top",
            fontsize=9.2,
            linespacing=1.35,
            bbox={"boxstyle": "round,pad=0.5", "facecolor": "#e6f7ff", "edgecolor": "#00a6d6"},
        )
        annotation_axis.text(
            0.345,
            0.98,
            outcomes,
            transform=annotation_axis.transAxes,
            va="top",
            fontsize=9.2,
            linespacing=1.25,
            bbox={"boxstyle": "round,pad=0.5", "facecolor": "#fff5f0", "edgecolor": "#d94801"},
        )
        annotation_axis.text(
            0.69,
            0.98,
            state,
            transform=annotation_axis.transAxes,
            va="top",
            fontsize=9.2,
            linespacing=1.22,
            bbox={"boxstyle": "round,pad=0.5", "facecolor": "#f7fcf5", "edgecolor": "#238b45"},
        )
        annotation_axis.legend(
            handles=[
                Patch(facecolor=color, edgecolor=color, alpha=0.35, label=label)
                for _, color, label in EVENT_OVERLAYS
            ],
            loc="lower center",
            bbox_to_anchor=(0.5, -0.02),
            ncol=5,
            fontsize=8.5,
            frameon=False,
        )
        return []

    animation = FuncAnimation(figure, update, frames=len(frames), blit=False)
    animation.save(
        OUTPUT_DIR / "connected_island_evolution_slow.gif",
        writer=PillowWriter(fps=GIF_FPS),
        dpi=105,
    )
    plt.close(figure)


def main() -> None:
    global ROAD_WIDTH_M, OUTPUT_DIR, GIF_FPS
    args = parse_args()
    if not np.isfinite(args.road_width) or args.road_width <= 0:
        raise ValueError("--road-width must be a finite positive number of meters.")
    if not np.isfinite(args.gif_fps) or args.gif_fps <= 0:
        raise ValueError("--gif-fps must be a finite positive number.")
    ROAD_WIDTH_M = float(args.road_width)
    OUTPUT_DIR = args.output_dir.resolve()
    GIF_FPS = float(args.gif_fps)

    input_dir, elevation_files, dune_files, input_validation = prepare_connected_inputs()
    write_manifest(input_dir, elevation_files, dune_files, input_validation, "running")
    cascade = None
    try:
        cascade, completed_year, stop_reason, event_log = run_model(
            input_dir, elevation_files, dune_files
        )
        metrics = extract_metrics(cascade, completed_year)
        np.savez_compressed(OUTPUT_DIR / "cascade_run.npz", cascade=cascade)
        np.savez_compressed(OUTPUT_DIR / "annual_metrics.npz", **metrics)
        write_event_log(event_log)
        write_domain_summary(metrics, cascade)
        report = validate_run(metrics, cascade)
        plot_connected_metrics(metrics)
        plot_connected_states(cascade, metrics)
        print("Rendering connected-island GIF...", flush=True)
        make_connected_gif(cascade, metrics)
        write_manifest(
            input_dir,
            elevation_files,
            dune_files,
            input_validation,
            "complete",
            int(metrics["year"][-1]),
            stop_reason,
        )
        print("\nConnected-island summary", flush=True)
        print(f"  final year: {metrics['year'][-1]}", flush=True)
        print(f"  active roads: {sum(not bool(v) for v in cascade.road_break)}/20", flush=True)
        print(f"  nourishment events: {int(np.sum(metrics['nourishment']))}", flush=True)
        print(f"  historical relocations: {int(np.sum(metrics['forced_relocation']))}", flush=True)
        print(f"  native relocations: {int(np.sum(metrics['triggered_relocation']))}", flush=True)
        print(f"  incomplete relocations: {int(np.sum(metrics['incomplete_relocation']))}", flush=True)
        print(f"  validation passed: {report['passed']}", flush=True)
        print(f"  outputs: {OUTPUT_DIR}", flush=True)
    except Exception as error:
        completed = None
        if cascade is not None:
            completed = min(len(b.x_s_TS) for b in cascade.barrier3d) - 1
        write_manifest(
            input_dir,
            elevation_files,
            dune_files,
            input_validation,
            "failed",
            completed,
            f"{type(error).__name__}: {error}",
        )
        raise


if __name__ == "__main__":
    main()
