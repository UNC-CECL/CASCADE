"""Create visual checks of the saved full-run nourishment comparison.

These tests read the existing BeachDuneManager and modified RoadwayManager NPZ
files. They do not run CASCADE and do not modify model or manager code. A passed
visualization test means the saved variables were valid and the diagnostic file
was created; the two full-manager histories are not required to be identical.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
import numpy as np
import pytest
from matplotlib.colors import BoundaryNorm, ListedColormap, TwoSlopeNorm
from matplotlib.patches import Patch

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


START_YEAR = 1992
END_YEAR = 2007
COMPARISON_DOMAINS = list(range(111, 120))
REAL_DOMAINS = list(range(80, 120))
FIRST_REAL_DOMAIN = 80
START_REAL_INDEX = 71

SOURCE_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = SOURCE_ROOT.parent
COMPARISON_ROOT = (
    WORKSPACE_ROOT
    / "CASCADE"
    / "output"
    / "roadway_nourishment"
    / "historical_manager_comparison"
)
BEACH_DUNE_NPZ = (
    COMPARISON_ROOT
    / "beach_dune_manager"
    / "PEA_1992_2007_HistoricalNourishment_BeachDuneManager_Hs2p0_Berm1p7"
    / "PEA_1992_2007_HistoricalNourishment_BeachDuneManager_Hs2p0_Berm1p7.npz"
)
ROADWAY_NPZ = (
    COMPARISON_ROOT
    / "roadway_manager"
    / "PEA_1992_2007_HistoricalNourishment_RoadwayManager_Hs2p0_Berm1p7"
    / "PEA_1992_2007_HistoricalNourishment_RoadwayManager_Hs2p0_Berm1p7.npz"
)
OUTPUT_DIR = COMPARISON_ROOT / "comparison_plots" / "tested_parameters"


def real_domain_to_index(domain: int) -> int:
    return START_REAL_INDEX + domain - FIRST_REAL_DOMAIN


def load_saved_cascade(path: Path):
    if not path.is_file():
        pytest.fail(f"Saved comparison run does not exist: {path}")
    with np.load(path, allow_pickle=True) as archive:
        return archive["cascade"].item()


@pytest.fixture(scope="module")
def saved_runs():
    return {
        "BeachDuneManager": load_saved_cascade(BEACH_DUNE_NPZ),
        "RoadwayManager": load_saved_cascade(ROADWAY_NPZ),
    }


def time_labels(number_of_states: int) -> list[str]:
    labels = ["Initial"] + [str(year) for year in range(START_YEAR, END_YEAR + 1)]
    assert len(labels) == number_of_states
    return labels


def barrier_matrix(cascade, attribute: str, scale: float = 1.0) -> np.ndarray:
    return np.vstack(
        [
            np.asarray(
                getattr(cascade.barrier3d[real_domain_to_index(domain)], attribute),
                dtype=float,
            )
            * scale
            for domain in COMPARISON_DOMAINS
        ]
    )


def manager_matrix(cascade, collection: str, attribute: str) -> np.ndarray:
    managers = getattr(cascade, collection)
    return np.vstack(
        [
            np.asarray(
                getattr(managers[real_domain_to_index(domain)], attribute),
                dtype=float,
            )
            for domain in COMPARISON_DOMAINS
        ]
    )


def comparison_data(saved_runs, source: str, attribute: str, scale: float = 1.0):
    beach = saved_runs["BeachDuneManager"]
    roadway = saved_runs["RoadwayManager"]
    if source == "barrier":
        beach_values = barrier_matrix(beach, attribute, scale)
        road_values = barrier_matrix(roadway, attribute, scale)
    elif source == "manager":
        beach_values = manager_matrix(beach, "nourishments", attribute) * scale
        road_values = manager_matrix(roadway, "roadways", attribute) * scale
    else:
        raise ValueError(f"Unknown data source: {source}")

    beach_events = manager_matrix(beach, "nourishments", "_nourishment_TS").astype(
        bool
    )
    road_events = manager_matrix(roadway, "roadways", "_nourishment_TS").astype(bool)
    return beach_values, road_values, beach_events, road_events


def validate_continuous_matrices(beach_values, road_values) -> None:
    expected_shape = (
        len(COMPARISON_DOMAINS),
        END_YEAR - START_YEAR + 2,
    )
    assert beach_values.shape == expected_shape
    assert road_values.shape == expected_shape
    assert np.isfinite(beach_values[:, 0]).all()
    assert np.isfinite(road_values[:, 0]).all()


def write_continuous_csv(
    path: Path,
    variable: str,
    beach_values: np.ndarray,
    road_values: np.ndarray,
    beach_events: np.ndarray,
    road_events: np.ndarray,
) -> None:
    labels = time_labels(beach_values.shape[1])
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=[
                "domain",
                "state_time",
                f"beach_dune_manager_{variable}",
                f"roadway_manager_{variable}",
                "roadway_minus_beach_dune",
                "beach_dune_manager_nourished",
                "roadway_manager_nourished",
            ],
        )
        writer.writeheader()
        for row, domain in enumerate(COMPARISON_DOMAINS):
            for column, label in enumerate(labels):
                beach_value = beach_values[row, column]
                road_value = road_values[row, column]
                writer.writerow(
                    {
                        "domain": domain,
                        "state_time": label,
                        f"beach_dune_manager_{variable}": beach_value,
                        f"roadway_manager_{variable}": road_value,
                        "roadway_minus_beach_dune": road_value - beach_value,
                        "beach_dune_manager_nourished": bool(
                            beach_events[row, column]
                        ),
                        "roadway_manager_nourished": bool(
                            road_events[row, column]
                        ),
                    }
                )


def plot_continuous_comparison(
    *,
    variable: str,
    display_name: str,
    units: str,
    beach_values: np.ndarray,
    road_values: np.ndarray,
    beach_events: np.ndarray,
    road_events: np.ndarray,
) -> tuple[Path, Path]:
    validate_continuous_matrices(beach_values, road_values)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    figure_path = OUTPUT_DIR / f"{variable}_beach_dune_vs_roadway.png"
    csv_path = OUTPUT_DIR / f"{variable}_beach_dune_vs_roadway.csv"
    write_continuous_csv(
        csv_path,
        variable,
        beach_values,
        road_values,
        beach_events,
        road_events,
    )

    combined = np.concatenate(
        [beach_values[np.isfinite(beach_values)], road_values[np.isfinite(road_values)]]
    )
    value_min = float(np.min(combined))
    value_max = float(np.max(combined))
    if np.isclose(value_min, value_max):
        value_max = value_min + 1.0

    difference = road_values - beach_values
    finite_difference = difference[np.isfinite(difference)]
    difference_limit = float(np.max(np.abs(finite_difference)))
    if np.isclose(difference_limit, 0.0):
        difference_limit = 1.0

    figure, axes = plt.subplots(3, 1, figsize=(17, 10), sharex=True)
    panels = [
        (beach_values, beach_events, "Original BeachDuneManager"),
        (road_values, road_events, "Modified RoadwayManager"),
    ]
    main_image = None
    for axis, (values, events, title) in zip(axes[:2], panels):
        color_map = plt.get_cmap("viridis").copy()
        color_map.set_bad("#d9d9d9")
        main_image = axis.imshow(
            np.ma.masked_invalid(values),
            aspect="auto",
            interpolation="none",
            cmap=color_map,
            vmin=value_min,
            vmax=value_max,
        )
        event_rows, event_columns = np.where(events)
        axis.scatter(
            event_columns,
            event_rows,
            marker="o",
            s=32,
            facecolor="#ffd92f",
            edgecolor="black",
            linewidth=0.6,
        )
        axis.set_title(title)
        axis.set_ylabel("Domain")
        axis.set_yticks(np.arange(len(COMPARISON_DOMAINS)), COMPARISON_DOMAINS)

    difference_map = plt.get_cmap("RdBu_r").copy()
    difference_map.set_bad("#d9d9d9")
    difference_image = axes[2].imshow(
        np.ma.masked_invalid(difference),
        aspect="auto",
        interpolation="none",
        cmap=difference_map,
        norm=TwoSlopeNorm(
            vmin=-difference_limit,
            vcenter=0.0,
            vmax=difference_limit,
        ),
    )
    axes[2].set_title("Modified RoadwayManager minus original BeachDuneManager")
    axes[2].set_ylabel("Domain")
    axes[2].set_yticks(np.arange(len(COMPARISON_DOMAINS)), COMPARISON_DOMAINS)

    labels = time_labels(beach_values.shape[1])
    axes[2].set_xticks(np.arange(len(labels)), labels, rotation=45, ha="right")
    axes[2].set_xlabel("Saved state (initial condition, then end of model year)")
    figure.colorbar(
        main_image,
        ax=axes[:2],
        location="right",
        shrink=0.82,
        label=f"{display_name} ({units})",
        pad=0.02,
    )
    figure.colorbar(
        difference_image,
        ax=axes[2],
        location="right",
        shrink=0.9,
        label=f"Difference ({units})",
        pad=0.02,
    )
    axes[0].legend(
        handles=[
            Patch(facecolor="#ffd92f", edgecolor="black", label="Applied nourishment"),
            Patch(facecolor="#d9d9d9", label="Unavailable"),
        ],
        loc="upper left",
        bbox_to_anchor=(1.07, 1.0),
    )
    figure.suptitle(
        f"{display_name}: full historical-nourishment runs, domains 111–119",
        fontsize=14,
    )
    figure.subplots_adjust(left=0.07, right=0.82, bottom=0.14, top=0.91, hspace=0.34)
    figure.savefig(figure_path, dpi=180, bbox_inches="tight")
    plt.close(figure)
    assert figure_path.is_file() and figure_path.stat().st_size > 0
    assert csv_path.is_file() and csv_path.stat().st_size > 0
    return figure_path, csv_path


def test_visualize_full_run_shoreline_position(saved_runs):
    values = comparison_data(saved_runs, "barrier", "x_s_TS", scale=10.0)
    plot_continuous_comparison(
        variable="shoreline_position_m",
        display_name="Shoreline position",
        units="m",
        beach_values=values[0],
        road_values=values[1],
        beach_events=values[2],
        road_events=values[3],
    )


def test_visualize_full_run_shoreface_slope(saved_runs):
    values = comparison_data(saved_runs, "barrier", "s_sf_TS")
    assert (values[0][np.isfinite(values[0])] > 0).all()
    assert (values[1][np.isfinite(values[1])] > 0).all()
    plot_continuous_comparison(
        variable="shoreface_slope",
        display_name="Shoreface slope",
        units="dimensionless",
        beach_values=values[0],
        road_values=values[1],
        beach_events=values[2],
        road_events=values[3],
    )


def test_visualize_full_run_beach_width(saved_runs):
    values = comparison_data(saved_runs, "manager", "_beach_width")
    plot_continuous_comparison(
        variable="beach_width_m",
        display_name="Managed beach width",
        units="m",
        beach_values=values[0],
        road_values=values[1],
        beach_events=values[2],
        road_events=values[3],
    )


def test_visualize_full_run_nourishment_flags_and_applied_volume(saved_runs):
    beach_volume, road_volume, beach_events, road_events = comparison_data(
        saved_runs,
        "manager",
        "_nourishment_volume_TS",
    )
    assert np.array_equal(beach_volume > 0, beach_events)
    assert np.array_equal(road_volume > 0, road_events)
    assert int(beach_events.sum()) == 30
    assert int(road_events.sum()) == 23

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    figure_path = OUTPUT_DIR / "nourishment_flags_and_applied_volume.png"
    csv_path = OUTPUT_DIR / "nourishment_flags_and_applied_volume.csv"
    labels = time_labels(beach_volume.shape[1])
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=[
                "domain",
                "state_time",
                "beach_dune_manager_event_flag",
                "roadway_manager_event_flag",
                "beach_dune_manager_applied_m3_per_m",
                "roadway_manager_applied_m3_per_m",
            ],
        )
        writer.writeheader()
        for row, domain in enumerate(COMPARISON_DOMAINS):
            for column, label in enumerate(labels):
                writer.writerow(
                    {
                        "domain": domain,
                        "state_time": label,
                        "beach_dune_manager_event_flag": bool(
                            beach_events[row, column]
                        ),
                        "roadway_manager_event_flag": bool(road_events[row, column]),
                        "beach_dune_manager_applied_m3_per_m": beach_volume[
                            row, column
                        ],
                        "roadway_manager_applied_m3_per_m": road_volume[row, column],
                    }
                )

    maximum_volume = max(float(beach_volume.max()), float(road_volume.max()))
    status = np.zeros(beach_events.shape, dtype=int)
    status[beach_events & road_events] = 1
    status[beach_events & ~road_events] = 2
    status[~beach_events & road_events] = 3
    status_map = ListedColormap(["#f7f7f7", "#41ab5d", "#ff8c00", "#984ea3"])
    status_norm = BoundaryNorm(np.arange(-0.5, 4.5, 1), status_map.N)

    figure, axes = plt.subplots(3, 1, figsize=(17, 10), sharex=True)
    volume_image = None
    for axis, volume, manager_name, event_count in [
        (axes[0], beach_volume, "Original BeachDuneManager", 30),
        (axes[1], road_volume, "Modified RoadwayManager", 23),
    ]:
        volume_image = axis.imshow(
            volume,
            aspect="auto",
            interpolation="none",
            cmap="YlGnBu",
            vmin=0,
            vmax=maximum_volume,
        )
        axis.set_title(f"{manager_name}: {event_count} applied events")
        axis.set_ylabel("Domain")
        axis.set_yticks(np.arange(len(COMPARISON_DOMAINS)), COMPARISON_DOMAINS)
    axes[2].imshow(
        status,
        aspect="auto",
        interpolation="none",
        cmap=status_map,
        norm=status_norm,
    )
    axes[2].set_title("Nourishment-event flag comparison")
    axes[2].set_ylabel("Domain")
    axes[2].set_yticks(np.arange(len(COMPARISON_DOMAINS)), COMPARISON_DOMAINS)
    axes[2].set_xticks(np.arange(len(labels)), labels, rotation=45, ha="right")
    axes[2].set_xlabel("Applied during model year")
    figure.colorbar(
        volume_image,
        ax=axes[:2],
        location="right",
        shrink=0.82,
        pad=0.02,
        label="Applied nourishment volume (m³/m)",
    )
    axes[2].legend(
        handles=[
            Patch(facecolor="#f7f7f7", edgecolor="black", label="Neither applied"),
            Patch(facecolor="#41ab5d", label="Both applied"),
            Patch(facecolor="#ff8c00", label="BeachDune only"),
            Patch(facecolor="#984ea3", label="Roadway only"),
        ],
        loc="upper left",
        bbox_to_anchor=(1.005, 1.0),
    )
    figure.suptitle(
        "Nourishment event flags and applied volumes, domains 111–119",
        fontsize=14,
    )
    figure.subplots_adjust(left=0.07, right=0.82, bottom=0.14, top=0.91, hspace=0.34)
    figure.savefig(figure_path, dpi=180, bbox_inches="tight")
    plt.close(figure)
    assert figure_path.is_file() and figure_path.stat().st_size > 0
    assert csv_path.is_file() and csv_path.stat().st_size > 0


RELOCATION_DIAGNOSTICS = [
    ("triggered_relocation_TS", "Triggered relocation"),
    ("relocation_incomplete_TS", "Incomplete relocation"),
    ("historical_relocation_requested_TS", "Historical request"),
    ("forced_relocation_TS", "Forced relocation completed"),
]


def relocation_matrix(cascade, attribute: str) -> np.ndarray:
    active = np.asarray(cascade.roadway_management_module, dtype=bool)
    rows = []
    for domain in REAL_DOMAINS:
        index = real_domain_to_index(domain)
        manager = cascade.roadways[index]
        values = np.asarray(getattr(manager, attribute), dtype=float)
        if not active[index]:
            values = np.full(values.shape, np.nan)
        rows.append(values)
    return np.vstack(rows)


def test_visualize_full_run_relocation_diagnostic_arrays(saved_runs):
    matrices = {
        scenario: {
            attribute: relocation_matrix(cascade, attribute)
            for attribute, _ in RELOCATION_DIAGNOSTICS
        }
        for scenario, cascade in saved_runs.items()
    }
    for scenario in matrices:
        assert matrices[scenario]["triggered_relocation_TS"][
            np.isfinite(matrices[scenario]["triggered_relocation_TS"])
        ].sum() == pytest.approx(3.0)
        for attribute in [
            "relocation_incomplete_TS",
            "historical_relocation_requested_TS",
            "forced_relocation_TS",
        ]:
            assert matrices[scenario][attribute][
                np.isfinite(matrices[scenario][attribute])
            ].sum() == pytest.approx(0.0)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    figure_path = OUTPUT_DIR / "relocation_diagnostic_arrays.png"
    csv_path = OUTPUT_DIR / "relocation_diagnostic_arrays.csv"
    labels = time_labels(next(iter(matrices["RoadwayManager"].values())).shape[1])
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=["scenario", "diagnostic", "domain", "state_time", "value"],
        )
        writer.writeheader()
        for scenario, diagnostic_matrices in matrices.items():
            for attribute, values in diagnostic_matrices.items():
                for row, domain in enumerate(REAL_DOMAINS):
                    for column, label in enumerate(labels):
                        value = values[row, column]
                        writer.writerow(
                            {
                                "scenario": scenario,
                                "diagnostic": attribute,
                                "domain": domain,
                                "state_time": label,
                                "value": "not applicable" if np.isnan(value) else int(value),
                            }
                        )

    event_map = ListedColormap(["#f7f7f7", "#e31a1c"])
    event_map.set_bad("#d9d9d9")
    event_norm = BoundaryNorm([-0.5, 0.5, 1.5], event_map.N)
    figure, axes = plt.subplots(4, 2, figsize=(19, 15), sharex=True, sharey=True)
    scenarios = ["BeachDuneManager", "RoadwayManager"]
    scenario_titles = {
        "BeachDuneManager": "BeachDuneManager scenario (road domains 82–110)",
        "RoadwayManager": "RoadwayManager scenario (road domains 82–119)",
    }
    displayed_domains = [80, 90, 100, 110, 119]
    displayed_positions = [REAL_DOMAINS.index(domain) for domain in displayed_domains]
    for row, (attribute, display_name) in enumerate(RELOCATION_DIAGNOSTICS):
        for column, scenario in enumerate(scenarios):
            values = matrices[scenario][attribute]
            count = int(np.nansum(values))
            axes[row, column].imshow(
                np.ma.masked_invalid(values),
                aspect="auto",
                interpolation="none",
                cmap=event_map,
                norm=event_norm,
            )
            event_rows, event_columns = np.where(values == 1)
            for event_row, event_column in zip(event_rows, event_columns):
                axes[row, column].text(
                    event_column,
                    event_row,
                    str(REAL_DOMAINS[event_row]),
                    ha="center",
                    va="center",
                    color="white",
                    fontsize=6,
                    fontweight="bold",
                )
            axes[row, column].set_title(
                f"{scenario_titles[scenario]}: {display_name} (count={count})"
            )
            axes[row, column].set_yticks(displayed_positions, displayed_domains)
            axes[row, column].set_ylabel("Domain")
    for axis in axes[-1, :]:
        axis.set_xticks(np.arange(len(labels)), labels, rotation=45, ha="right")
        axis.set_xlabel("Saved state")
    axes[0, 1].legend(
        handles=[
            Patch(facecolor="#f7f7f7", edgecolor="black", label="No event"),
            Patch(facecolor="#e31a1c", label="Event"),
            Patch(facecolor="#d9d9d9", label="No roadway manager"),
        ],
        loc="upper left",
        bbox_to_anchor=(1.005, 1.0),
    )
    figure.suptitle(
        "Roadway relocation diagnostic arrays in both full-run scenarios\n"
        "All real domains 80–119",
        fontsize=15,
    )
    figure.subplots_adjust(left=0.06, right=0.87, bottom=0.09, top=0.93, hspace=0.3)
    figure.savefig(figure_path, dpi=180, bbox_inches="tight")
    plt.close(figure)
    assert figure_path.is_file() and figure_path.stat().st_size > 0
    assert csv_path.is_file() and csv_path.stat().st_size > 0
