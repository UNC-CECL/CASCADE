#!/usr/bin/env python3
"""Run and plot an unprescribed synthetic RoadwayManager comparison.

Both cases use the same real Barrier3D input files and normal ``Cascade.update``
loop. The only difference is whether the new road-setback dune-growth rule is
enabled. No nourishment, historical relocation, storm suppression, shoreline
movement, dune movement, or time-varying road state is prescribed.
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
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import Rectangle

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SOURCE_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = SOURCE_ROOT.parent
sys.path.insert(0, str(SOURCE_ROOT))

from cascade import Cascade  # noqa: E402


INPUT_ROOT = SOURCE_ROOT / "tests" / "test_human_dynamics"
INPUT_FILES = (
    "StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
    "b3d_pt75_3284yrs_low-elevations.csv",
    "pathways-dunes.npy",
    "roadway-parameters.yaml",
)
DEFAULT_GROWTH_PARAM_SOURCE = INPUT_ROOT / "growthparam_1000dam.npy"
GROWTH_PARAM_SOURCE = DEFAULT_GROWTH_PARAM_SOURCE
TIME_STEP_COUNT = 180
INITIAL_BEACH_WIDTH_M = 30.0

OUTPUT_DIR = (
    WORKSPACE_ROOT
    / "CASCADE"
    / "output"
    / "roadway_dune_growth_management"
    / "synthetic_actual_barrier3d_natural_evolution"
)
PNG_PATH = OUTPUT_DIR / "controlled_style_growth_management_with_dmax.png"
GIF_PATH = OUTPUT_DIR / "controlled_style_growth_management_with_dmax_slow.gif"
CSV_PATH = OUTPUT_DIR / "growth_management_comparison.csv"
MANIFEST_PATH = OUTPUT_DIR / "growth_management_comparison_manifest.json"

SCENARIOS = {
    "Original growth logic (rule OFF)": False,
    "New setback rule (rule ON)": True,
}
COLORS = {
    "Original growth logic (rule OFF)": "#238b45",
    "New setback rule (rule ON)": "#d95f0e",
}

COMMON_PARAMETERS = {
    "storm_file": "StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
    "elevation_file": "b3d_pt75_3284yrs_low-elevations.csv",
    "dune_file": "pathways-dunes.npy",
    "parameter_file": "roadway-parameters.yaml",
    "wave_height": 1,
    "wave_period": 7,
    "wave_asymmetry": 0.8,
    "wave_angle_high_fraction": 0.2,
    "sea_level_rise_rate": 0.004,
    "sea_level_rise_constant": True,
    "background_erosion": 0.0,
    "alongshore_section_count": 1,
    "time_step_count": TIME_STEP_COUNT,
    "min_dune_growth_rate": 0.55,
    "max_dune_growth_rate": 0.95,
    "num_cores": 1,
    "roadway_management_module": True,
    "alongshore_transport_module": False,
    "beach_nourishment_module": False,
    "community_economics_module": False,
    "road_ele": 1.2,
    "road_width": 20,
    "road_setback": 20,
    "dune_design_elevation": 3.2,
    "dune_minimum_elevation": 1.7,
    "nourishment_interval": None,
    "outwash_module": False,
}


def copy_inputs(run_directory: Path) -> None:
    """Copy inputs and select the requested growth-parameter file."""

    run_directory.mkdir(parents=True, exist_ok=True)
    for filename in INPUT_FILES:
        shutil.copy2(INPUT_ROOT / filename, run_directory / filename)
    shutil.copy2(GROWTH_PARAM_SOURCE, run_directory / GROWTH_PARAM_SOURCE.name)

    parameter_path = run_directory / "roadway-parameters.yaml"
    lines = parameter_path.read_text().splitlines()
    matches = [
        index for index, line in enumerate(lines) if line.startswith("growth_param_file:")
    ]
    if len(matches) != 1:
        raise RuntimeError(
            f"Expected one growth_param_file entry in {parameter_path}; found {matches}"
        )
    lines[matches[0]] = f"growth_param_file: {GROWTH_PARAM_SOURCE.name}"
    parameter_path.write_text("\n".join(lines) + "\n")


def load_growth_parameter_input(cell_count: int) -> np.ndarray:
    """Load exactly one initial alongshore growth array from the requested file."""

    if GROWTH_PARAM_SOURCE.suffix == ".npy":
        values = np.load(GROWTH_PARAM_SOURCE)
    else:
        values = np.loadtxt(GROWTH_PARAM_SOURCE)
    values = np.asarray(values, dtype=float).reshape(-1)
    if values.size < cell_count:
        raise ValueError(
            f"{GROWTH_PARAM_SOURCE} contains {values.size} values; "
            f"the Barrier3D segment requires {cell_count}."
        )
    return values[:cell_count].reshape(1, cell_count)


def run_scenario(label: str, enabled: bool) -> Cascade:
    """Run one case with no state changes prescribed after initialization."""

    run_directory = OUTPUT_DIR / "run_inputs" / ("rule_on" if enabled else "rule_off")
    copy_inputs(run_directory)
    cascade = Cascade(
        str(run_directory),
        name=(
            "synthetic_road_dune_growth_rule_on"
            if enabled
            else "synthetic_original_growth_logic"
        ),
        road_dune_growth_management_on=enabled,
        **COMMON_PARAMETERS,
    )
    # CASCADE's initialize_equal() always sets GrowthParamStart=False and creates
    # random growth parameters from rmin/rmax. Apply the explicitly requested file
    # once here, before the first update, so Barrier3D and RoadwayManager share the
    # exact same initial and restoration array. Nothing is assigned in the loop.
    barrier = cascade.barrier3d[0]
    roadway = cascade.roadways[0]
    requested_growth_param = load_growth_parameter_input(barrier.growthparam.shape[1])
    barrier.growthparam = requested_growth_param.copy()
    roadway._original_growth_param = requested_growth_param.copy()
    roadway._growth_params[0] = requested_growth_param.copy()

    for _ in range(TIME_STEP_COUNT - 1):
        cascade.update()
        if cascade.b3d_break:
            break

    np.savez(
        OUTPUT_DIR / ("rule_on.npz" if enabled else "rule_off.npz"),
        cascade=[cascade],
    )
    print(f"{label}: completed {len(cascade.barrier3d[0].x_s_TS)} states")
    return cascade


def growth_series(roadway, count: int) -> tuple[np.ndarray, np.ndarray]:
    means = np.full(count, np.nan)
    zero_fraction = np.full(count, np.nan)
    for index, values in enumerate(roadway._growth_params[:count]):
        if isinstance(values, np.ndarray):
            values = np.asarray(values, dtype=float)
            means[index] = np.mean(values)
            zero_fraction[index] = np.mean(values == 0.0)
    return means, zero_fraction


def extract(cascade: Cascade) -> dict[str, np.ndarray]:
    barrier = cascade.barrier3d[0]
    roadway = cascade.roadways[0]
    count = len(barrier.x_s_TS)
    shoreline_change = np.asarray(barrier.ShorelineChangeTS[:count], dtype=float)
    shoreline_m = np.asarray(barrier.x_s_TS[:count], dtype=float) * 10.0
    initial_dune_grid_m = (
        np.floor(barrier.x_s_TS[0] + INITIAL_BEACH_WIDTH_M / 10.0) * 10.0
    )
    dune_grid_m = initial_dune_grid_m + np.cumsum(-shoreline_change) * 10.0
    road_setback_m = np.asarray(roadway._road_setback_TS[:count], dtype=float)
    road_width_m = np.asarray(roadway._road_width_TS[:count], dtype=float)
    road_position_m = dune_grid_m + road_setback_m
    road_position_m[road_width_m <= 0] = np.nan
    growth_mean, growth_zero_fraction = growth_series(roadway, count)

    dune_crest_mean_m = np.full(count, np.nan)
    front_row_min_m = np.full(count, np.nan)
    front_row_mean_m = np.full(count, np.nan)
    front_row_max_m = np.full(count, np.nan)
    fraction_front_row_above_dmax = np.full(count, np.nan)
    dmax_m = float(barrier.Dmax) * 10.0
    for index in range(count):
        dune_height_m = np.asarray(barrier.DuneDomain[index], dtype=float) * 10.0
        dune_crest_mean_m[index] = np.mean(np.max(dune_height_m, axis=1)) + (
            barrier.BermEl * 10.0
        )
        front_row = dune_height_m[:, 0]
        front_row_min_m[index] = np.min(front_row)
        front_row_mean_m[index] = np.mean(front_row)
        front_row_max_m[index] = np.max(front_row)
        fraction_front_row_above_dmax[index] = np.mean(front_row > dmax_m)

    return {
        "time": np.arange(count),
        "shoreline_m": shoreline_m,
        "dune_grid_m": dune_grid_m,
        "road_position_m": road_position_m,
        "road_width_m": road_width_m,
        "road_setback_m": road_setback_m,
        "growth_mean": growth_mean,
        "growth_zero_fraction": growth_zero_fraction,
        "dune_crest_mean_m": dune_crest_mean_m,
        "front_row_min_m": front_row_min_m,
        "front_row_mean_m": front_row_mean_m,
        "front_row_max_m": front_row_max_m,
        "fraction_front_row_above_dmax": fraction_front_row_above_dmax,
        "dmax_m": np.full(count, dmax_m),
        "management_active": np.asarray(
            roadway.road_dune_growth_management_TS[:count], dtype=bool
        ),
        "triggered_relocation": np.asarray(
            roadway.triggered_relocation_TS[:count], dtype=bool
        ),
        "incomplete_relocation": np.asarray(
            roadway.relocation_incomplete_TS[:count], dtype=bool
        ),
        "historical_request": np.asarray(
            roadway.historical_relocation_requested_TS[:count], dtype=bool
        ),
        "forced_relocation": np.asarray(
            roadway.forced_relocation_TS[:count], dtype=bool
        ),
        "nourishment": np.asarray(roadway.nourishment_TS[:count], dtype=bool),
        "dunes_rebuilt": np.asarray(roadway._dunes_rebuilt_TS[:count], dtype=bool),
        "storm_count": np.asarray(barrier._StormCount[:count], dtype=int),
        "qow": np.asarray(barrier.QowTS[:count], dtype=float),
        "average_barrier_width_m": np.asarray(
            barrier.InteriorWidth_AvgTS[:count], dtype=float
        )
        * 10.0,
    }


def value_at(values: dict[str, np.ndarray], key: str, index: int):
    return values[key][min(index, len(values[key]) - 1)]


def draw_plan(axis, label: str, values: dict[str, np.ndarray], index: int, limits) -> None:
    actual_index = min(index, len(values["time"]) - 1)
    shoreline = float(values["shoreline_m"][actual_index])
    dune = float(values["dune_grid_m"][actual_index])
    road = float(values["road_position_m"][actual_index])
    road_width = float(values["road_width_m"][actual_index])
    barrier_width = float(values["average_barrier_width_m"][actual_index])
    bay = dune + barrier_width
    x_min, x_max = limits

    axis.axvspan(x_min, shoreline, color="#9ecae1", alpha=0.95, label="Ocean")
    axis.axvspan(shoreline, dune, color="#fdd49e", alpha=0.95, label="Beach")
    axis.axvspan(dune, bay, color="#c7e9c0", alpha=0.95, label="Barrier interior")
    axis.axvspan(bay, x_max, color="#a6bddb", alpha=0.75, label="Back-barrier water")
    axis.axvline(shoreline, color="#08519c", linewidth=2.5, label="Shoreline")
    axis.axvspan(dune - 1.2, dune + 1.2, color="#8c510a", label="Dune grid")
    if np.isfinite(road) and road_width > 0:
        axis.add_patch(
            Rectangle(
                (road, 17),
                road_width,
                16,
                facecolor="#636363",
                edgecolor="black",
                label="Road",
            )
        )

    active = bool(values["management_active"][actual_index])
    growth = float(values["growth_mean"][actual_index])
    zero_fraction = float(values["growth_zero_fraction"][actual_index])
    status = "ACTIVE: growthparam = 0" if active else "inactive: original logic"
    if actual_index < index:
        status += f"\nrun ended at index {actual_index}"
    axis.text(
        0.02,
        0.97,
        f"Setback: {values['road_setback_m'][actual_index]:.1f} m\n"
        f"Rule: {status}\n"
        f"Mean growthparam: {growth:.3f}\n"
        f"Zero-valued cells: {zero_fraction * 100:.1f}%\n"
        f"Front-row dunes: {values['front_row_min_m'][actual_index]:.2f}–"
        f"{values['front_row_max_m'][actual_index]:.2f} m above berm\n"
        f"Dmax: {values['dmax_m'][actual_index]:.2f} m above berm\n"
        f"Storms: {values['storm_count'][actual_index]} | Qow: {values['qow'][actual_index]:.1f}",
        transform=axis.transAxes,
        va="top",
        fontsize=8.5,
        bbox={
            "boxstyle": "round,pad=0.4",
            "facecolor": "#fff7bc" if active else "white",
            "edgecolor": "#d95f0e" if active else "#238b45",
            "linewidth": 2,
            "alpha": 0.96,
        },
    )
    axis.set_xlim(x_min, x_max)
    axis.set_ylim(0, 50)
    axis.set_yticks([])
    axis.set_xlabel("Cross-shore position (m; landward →)")
    axis.set_title(label, fontweight="bold")


def make_figure(data: dict[str, dict[str, np.ndarray]], index: int):
    figure, axes = plt.subplots(3, 2, figsize=(15, 12))
    all_shorelines = np.concatenate([values["shoreline_m"] for values in data.values()])
    all_dunes = np.concatenate([values["dune_grid_m"] for values in data.values()])
    all_bays = np.concatenate(
        [
            values["dune_grid_m"] + values["average_barrier_width_m"]
            for values in data.values()
        ]
    )
    limits = (
        float(all_shorelines.min() - 20),
        float(max(all_dunes.max() + 60, all_bays.max() + 15)),
    )

    labels = list(data)
    draw_plan(axes[0, 0], labels[0], data[labels[0]], index, limits)
    draw_plan(axes[0, 1], labels[1], data[labels[1]], index, limits)
    axes[0, 1].legend(loc="lower right", fontsize=7, ncol=2)

    for label, values in data.items():
        color = COLORS[label]
        stop = min(index + 1, len(values["time"]))
        time = values["time"][:stop]
        axes[1, 0].plot(
            time,
            values["road_setback_m"][:stop],
            color=color,
            linewidth=2,
            label=label,
        )
        relocation = np.flatnonzero(values["triggered_relocation"][:stop])
        axes[1, 0].scatter(
            relocation,
            values["road_setback_m"][relocation],
            marker="*",
            s=110,
            color=color,
            edgecolor="black",
            zorder=5,
        )
        active = np.flatnonzero(values["management_active"][:stop])
        axes[1, 0].scatter(
            active,
            values["road_setback_m"][active],
            s=16,
            color="#e31a1c",
            zorder=4,
        )

        axes[1, 1].plot(
            time,
            values["growth_mean"][:stop],
            color=color,
            linewidth=2.2,
            label=f"{label}: mean growthparam",
        )
        axes[2, 0].fill_between(
            time,
            values["front_row_min_m"][:stop],
            values["front_row_max_m"][:stop],
            color=color,
            alpha=0.14,
        )
        axes[2, 0].plot(
            time,
            values["front_row_mean_m"][:stop],
            color=color,
            linewidth=2,
            label=f"{label}: mean front-row height",
        )
        rebuilt = np.flatnonzero(values["dunes_rebuilt"][:stop])
        axes[2, 0].scatter(
            rebuilt,
            values["front_row_mean_m"][rebuilt],
            marker="P",
            s=65,
            color=color,
            edgecolor="black",
            zorder=5,
        )

        axes[2, 1].plot(
            time,
            values["fraction_front_row_above_dmax"][:stop],
            color=color,
            linestyle="--",
            linewidth=2,
            label=f"{label}: front-row cells > Dmax",
        )
        axes[2, 1].plot(
            time,
            values["growth_zero_fraction"][:stop],
            color=color,
            linestyle=":",
            linewidth=1.8,
            label=f"{label}: growthparam = 0",
        )

    max_count = max(len(values["time"]) for values in data.values())
    axes[1, 0].axhline(
        20.0,
        color="black",
        linestyle="--",
        linewidth=1.3,
        label="20 m threshold",
    )
    axes[1, 0].axhspan(0, 20, color="#fee8c8", alpha=0.35)
    axes[1, 0].set_xlim(0, max_count - 1)
    axes[1, 0].set_xlabel("Actual Barrier3D time index")
    axes[1, 0].set_ylabel("Road setback (m)")
    axes[1, 0].set_title("Calculated road setback (red dots = new rule active)")
    axes[1, 0].grid(alpha=0.25)
    axes[1, 0].legend(fontsize=7)

    axes[1, 1].set_xlim(0, max_count - 1)
    axes[1, 1].set_ylim(-0.03, 1.03)
    axes[1, 1].set_xlabel("Actual Barrier3D time index")
    axes[1, 1].set_ylabel("Mean growth parameter")
    axes[1, 1].set_title("Barrier3D mean dune growth parameter")
    axes[1, 1].grid(alpha=0.25)
    axes[1, 1].legend(fontsize=7.2)

    dmax_m = float(next(iter(data.values()))["dmax_m"][0])
    axes[2, 0].axhline(
        dmax_m,
        color="black",
        linestyle="--",
        linewidth=1.3,
        label=f"Dmax = {dmax_m:.2f} m above berm",
    )
    axes[2, 0].set_xlim(0, max_count - 1)
    axes[2, 0].set_xlabel("Actual Barrier3D time index")
    axes[2, 0].set_ylabel("Front-row dune height above berm (m)")
    axes[2, 0].set_title(
        "Exact height used by original growthparam logic\n"
        "(shading = alongshore min–max; P = dune rebuild)"
    )
    axes[2, 0].grid(alpha=0.25)
    axes[2, 0].legend(fontsize=6.8)

    axes[2, 1].set_xlim(0, max_count - 1)
    axes[2, 1].set_ylim(-0.03, 1.03)
    axes[2, 1].set_xlabel("Actual Barrier3D time index")
    axes[2, 1].set_ylabel("Fraction of alongshore cells")
    axes[2, 1].set_title(
        "Why growthparam is zero\n"
        "(original Dmax condition versus applied zero values)"
    )
    axes[2, 1].grid(alpha=0.25)
    axes[2, 1].legend(fontsize=6.4)

    event_messages = []
    for label, values in data.items():
        actual_index = min(index, len(values["time"]) - 1)
        if values["management_active"][actual_index]:
            event_messages.append(f"{label}: setback < 20 m; growthparam zero")
        if values["triggered_relocation"][actual_index]:
            event_messages.append(f"{label}: model-triggered road relocation")
        if values["dunes_rebuilt"][actual_index]:
            event_messages.append(f"{label}: original roadway dune rebuild")
        if values["storm_count"][actual_index] > 0:
            event_messages.append(
                f"{label}: {values['storm_count'][actual_index]} calculated storm(s)"
            )
    subtitle = (
        " | ".join(event_messages)
        if event_messages
        else "All evolution calculated by Barrier3D and RoadwayManager"
    )
    figure.suptitle(
        f"Synthetic actual-Barrier3D dune-growth comparison — time index {index}\n{subtitle}",
        fontsize=13,
        fontweight="bold",
    )
    figure.tight_layout(rect=(0, 0, 1, 0.92))
    return figure


def save_csv(data: dict[str, dict[str, np.ndarray]]) -> None:
    fields = (
        "time",
        "road_setback_m",
        "growth_mean",
        "growth_zero_fraction",
        "dune_crest_mean_m",
        "front_row_min_m",
        "front_row_mean_m",
        "front_row_max_m",
        "fraction_front_row_above_dmax",
        "dmax_m",
        "management_active",
        "triggered_relocation",
        "dunes_rebuilt",
        "storm_count",
        "qow",
    )
    with CSV_PATH.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("scenario", *fields))
        for label, values in data.items():
            for index in range(len(values["time"])):
                writer.writerow((label, *(values[field][index] for field in fields)))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_unprescribed(data: dict[str, dict[str, np.ndarray]]) -> None:
    for label, values in data.items():
        if values["nourishment"].any():
            raise RuntimeError(f"Unexpected nourishment in {label}")
        if values["historical_request"].any() or values["forced_relocation"].any():
            raise RuntimeError(f"Unexpected requested relocation in {label}")


def configure_paths(output_directory: Path, growth_param_source: Path) -> None:
    global OUTPUT_DIR
    global PNG_PATH
    global GIF_PATH
    global CSV_PATH
    global MANIFEST_PATH
    global GROWTH_PARAM_SOURCE

    OUTPUT_DIR = output_directory.expanduser().resolve()
    GROWTH_PARAM_SOURCE = growth_param_source.expanduser().resolve()
    PNG_PATH = OUTPUT_DIR / "controlled_style_growth_management_with_dmax.png"
    GIF_PATH = OUTPUT_DIR / "controlled_style_growth_management_with_dmax_slow.gif"
    CSV_PATH = OUTPUT_DIR / "growth_management_comparison.csv"
    MANIFEST_PATH = OUTPUT_DIR / "growth_management_comparison_manifest.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--growthparam",
        type=Path,
        default=DEFAULT_GROWTH_PARAM_SOURCE,
        help="Barrier3D growth-parameter input file.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=OUTPUT_DIR,
        help="Separate directory for saved runs and plots.",
    )
    parser.add_argument(
        "--reuse-saved",
        action="store_true",
        help="Load rule_off.npz and rule_on.npz instead of rerunning the model.",
    )
    return parser.parse_args()


def main() -> None:
    arguments = parse_args()
    configure_paths(arguments.output_dir, arguments.growthparam)
    if not GROWTH_PARAM_SOURCE.is_file():
        raise FileNotFoundError(GROWTH_PARAM_SOURCE)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    if arguments.reuse_saved:
        runs = {}
        for label, enabled in SCENARIOS.items():
            path = OUTPUT_DIR / ("rule_on.npz" if enabled else "rule_off.npz")
            if not path.is_file():
                raise FileNotFoundError(path)
            with np.load(path, allow_pickle=True) as archive:
                runs[label] = archive["cascade"].item()
    else:
        runs = {
            label: run_scenario(label, enabled)
            for label, enabled in SCENARIOS.items()
        }
    data = {label: extract(cascade) for label, cascade in runs.items()}
    validate_unprescribed(data)
    save_csv(data)

    max_count = max(len(values["time"]) for values in data.values())
    final_figure = make_figure(data, max_count - 1)
    final_figure.savefig(PNG_PATH, dpi=190, bbox_inches="tight")
    plt.close(final_figure)

    animation_figure = plt.figure(figsize=(15, 12))

    def draw_fixed(frame_index: int):
        animation_figure.clear()
        frame = make_figure(data, frame_index)
        frame.canvas.draw()
        image = np.asarray(frame.canvas.buffer_rgba()).copy()
        plt.close(frame)
        axis = animation_figure.add_axes((0, 0, 1, 1))
        axis.imshow(image)
        axis.axis("off")
        return (axis,)

    frame_indices = []
    event_indices = set()
    for values in data.values():
        event_indices.update(np.flatnonzero(values["management_active"]).tolist())
        event_indices.update(np.flatnonzero(values["triggered_relocation"]).tolist())
    for index in range(max_count):
        frame_indices.extend([index] * (2 if index in event_indices else 1))
    animation = FuncAnimation(
        animation_figure,
        draw_fixed,
        frames=frame_indices,
        interval=1000,
        blit=False,
    )
    animation.save(GIF_PATH, writer=PillowWriter(fps=1), dpi=95)
    plt.close(animation_figure)

    manifest = {
        "experiment": "synthetic actual-Barrier3D natural evolution",
        "source_root": str(SOURCE_ROOT),
        "output_directory": str(OUTPUT_DIR),
        "time_step_count_requested": TIME_STEP_COUNT,
        "only_scenario_difference": "road_dune_growth_management_on",
        "scenario_values": SCENARIOS,
        "common_parameters": COMMON_PARAMETERS,
        "prescribed_during_time_loop": [],
        "initial_growth_parameter_assignment": (
            "The requested file was loaded once immediately after Cascade "
            "construction because initialize_equal() forces GrowthParamStart=False."
        ),
        "nourishment_requests": 0,
        "historical_relocation_requests": 0,
        "storm_inputs_modified": False,
        "simulation_reused_for_plotting": arguments.reuse_saved,
        "input_sha256": {
            filename: sha256(INPUT_ROOT / filename) for filename in INPUT_FILES
        },
        "growth_parameter_input": {
            "path": str(GROWTH_PARAM_SOURCE),
            "copied_filename": GROWTH_PARAM_SOURCE.name,
            "sha256": sha256(GROWTH_PARAM_SOURCE),
        },
        "outputs": {
            "png": str(PNG_PATH),
            "gif": str(GIF_PATH),
            "csv": str(CSV_PATH),
            "rule_off_npz": str(OUTPUT_DIR / "rule_off.npz"),
            "rule_on_npz": str(OUTPUT_DIR / "rule_on.npz"),
        },
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n")
    print(PNG_PATH)
    print(GIF_PATH)
    print(CSV_PATH)
    print(MANIFEST_PATH)


if __name__ == "__main__":
    main()
