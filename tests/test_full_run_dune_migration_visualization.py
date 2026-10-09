"""Validate and visualize dune-migration state in the two saved full runs.

This test reads existing simulation output only. It does not rerun CASCADE or
change either management module. The figure distinguishes an active dune-
migration state from a dune line held in place by beach-width management, and
uses gray for years after a manager stopped recording state.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
import numpy as np
import pytest
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


START_YEAR = 1992
END_YEAR = 2007
DOMAINS = list(range(111, 120))
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
FIGURE_DIR = COMPARISON_ROOT / "comparison_plots" / "dune_migration"
FIGURE_PATH = FIGURE_DIR / "dune_migration_state_beach_dune_vs_roadway.png"
CSV_PATH = FIGURE_DIR / "dune_migration_state_beach_dune_vs_roadway.csv"


def real_domain_to_index(domain: int) -> int:
    return START_REAL_INDEX + domain - FIRST_REAL_DOMAIN


def load_saved_cascade(path: Path):
    if not path.is_file():
        pytest.fail(f"Saved comparison run does not exist: {path}")
    with np.load(path, allow_pickle=True) as archive:
        return archive["cascade"].item()


def manager_matrices(cascade, manager_attribute: str) -> tuple[np.ndarray, np.ndarray]:
    managers = getattr(cascade, manager_attribute)
    state_rows = []
    nourishment_rows = []
    for domain in DOMAINS:
        manager = managers[real_domain_to_index(domain)]
        state_rows.append(np.asarray(manager._dune_migration_on, dtype=float))
        nourishment_rows.append(np.asarray(manager._nourishment_TS, dtype=bool))
    return np.vstack(state_rows), np.vstack(nourishment_rows)


@pytest.fixture(scope="module")
def full_run_states():
    beach_dune = load_saved_cascade(BEACH_DUNE_NPZ)
    roadway = load_saved_cascade(ROADWAY_NPZ)
    beach_state, beach_nourishment = manager_matrices(
        beach_dune, "nourishments"
    )
    road_state, road_nourishment = manager_matrices(roadway, "roadways")
    return {
        "BeachDuneManager": (beach_state, beach_nourishment),
        "RoadwayManager": (road_state, road_nourishment),
    }


@pytest.mark.parametrize("manager_name", ["BeachDuneManager", "RoadwayManager"])
def test_full_run_dune_migration_states_are_boolean_or_unavailable(
    full_run_states,
    manager_name,
):
    """Every recorded state must be OFF (0), ON (1), or unavailable (NaN)."""

    states, _ = full_run_states[manager_name]
    assert states.shape == (len(DOMAINS), END_YEAR - START_YEAR + 2)
    assert np.isin(states[np.isfinite(states)], [0.0, 1.0]).all()


@pytest.mark.parametrize("manager_name", ["BeachDuneManager", "RoadwayManager"])
def test_applied_nourishment_leaves_dune_migration_off(
    full_run_states,
    manager_name,
):
    """At every applied historical event, dune migration must be OFF afterward."""

    states, nourishment = full_run_states[manager_name]
    assert nourishment.any(), f"No applied nourishment found for {manager_name}"
    np.testing.assert_array_equal(states[nourishment], 0.0)


def time_labels(number_of_states: int) -> list[str]:
    labels = ["Initial"]
    labels.extend(f"After {year}" for year in range(START_YEAR, END_YEAR + 1))
    assert len(labels) == number_of_states
    return labels


def write_state_csv(
    beach_state: np.ndarray,
    road_state: np.ndarray,
    beach_nourishment: np.ndarray,
    road_nourishment: np.ndarray,
) -> None:
    labels = time_labels(beach_state.shape[1])
    with CSV_PATH.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=[
                "domain",
                "state_time",
                "beach_dune_manager_dune_migration_on",
                "roadway_manager_dune_migration_on",
                "states_match",
                "beach_dune_manager_nourished",
                "roadway_manager_nourished",
            ],
        )
        writer.writeheader()
        for row, domain in enumerate(DOMAINS):
            for column, label in enumerate(labels):
                beach_value = beach_state[row, column]
                road_value = road_state[row, column]
                writer.writerow(
                    {
                        "domain": domain,
                        "state_time": label,
                        "beach_dune_manager_dune_migration_on": (
                            "unavailable" if np.isnan(beach_value) else int(beach_value)
                        ),
                        "roadway_manager_dune_migration_on": (
                            "unavailable" if np.isnan(road_value) else int(road_value)
                        ),
                        "states_match": (
                            "unavailable"
                            if np.isnan(beach_value) or np.isnan(road_value)
                            else bool(beach_value == road_value)
                        ),
                        "beach_dune_manager_nourished": bool(
                            beach_nourishment[row, column]
                        ),
                        "roadway_manager_nourished": bool(
                            road_nourishment[row, column]
                        ),
                    }
                )


def test_create_full_run_dune_migration_comparison_figure(full_run_states):
    """Create heatmaps from the exact arrays validated by the tests above."""

    beach_state, beach_nourishment = full_run_states["BeachDuneManager"]
    road_state, road_nourishment = full_run_states["RoadwayManager"]
    assert beach_state.shape == road_state.shape

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    write_state_csv(
        beach_state,
        road_state,
        beach_nourishment,
        road_nourishment,
    )

    state_cmap = ListedColormap(["#225ea8", "#41ab5d"])
    state_cmap.set_bad("#d9d9d9")
    state_norm = BoundaryNorm([-0.5, 0.5, 1.5], state_cmap.N)

    comparison = np.full(beach_state.shape, 4, dtype=int)
    available = np.isfinite(beach_state) & np.isfinite(road_state)
    comparison[available & (beach_state == 0) & (road_state == 0)] = 0
    comparison[available & (beach_state == 1) & (road_state == 1)] = 1
    comparison[available & (beach_state == 1) & (road_state == 0)] = 2
    comparison[available & (beach_state == 0) & (road_state == 1)] = 3
    comparison_cmap = ListedColormap(
        ["#225ea8", "#41ab5d", "#984ea3", "#ff8c00", "#d9d9d9"]
    )
    comparison_norm = BoundaryNorm(np.arange(-0.5, 5.5, 1), comparison_cmap.N)

    figure, axes = plt.subplots(3, 1, figsize=(17, 10), sharex=True)
    panels = [
        (beach_state, beach_nourishment, "Original BeachDuneManager"),
        (road_state, road_nourishment, "Modified RoadwayManager"),
    ]
    for axis, (states, nourishment, title) in zip(axes[:2], panels):
        axis.imshow(
            np.ma.masked_invalid(states),
            aspect="auto",
            interpolation="none",
            cmap=state_cmap,
            norm=state_norm,
        )
        event_rows, event_columns = np.where(nourishment)
        axis.scatter(
            event_columns,
            event_rows,
            marker="o",
            s=34,
            facecolor="#ffd92f",
            edgecolor="black",
            linewidth=0.6,
            label="Applied nourishment",
        )
        axis.set_title(title)
        axis.set_ylabel("Domain")
        axis.set_yticks(np.arange(len(DOMAINS)), labels=DOMAINS)
        axis.legend(loc="upper left", bbox_to_anchor=(1.005, 1.0))

    axes[2].imshow(
        comparison,
        aspect="auto",
        interpolation="none",
        cmap=comparison_cmap,
        norm=comparison_norm,
    )
    axes[2].set_title("State comparison")
    axes[2].set_ylabel("Domain")
    axes[2].set_yticks(np.arange(len(DOMAINS)), labels=DOMAINS)

    labels = time_labels(beach_state.shape[1])
    axes[2].set_xticks(np.arange(len(labels)), labels=labels, rotation=50, ha="right")
    axes[2].set_xlabel("Saved state (initial condition, then end of model year)")

    axes[0].legend(
        handles=[
            Patch(facecolor="#225ea8", label="Dune migration OFF"),
            Patch(facecolor="#41ab5d", label="Dune migration ON"),
            Patch(facecolor="#d9d9d9", label="Unavailable after management stops"),
            axes[0].collections[0],
        ],
        labels=[
            "Dune migration OFF",
            "Dune migration ON",
            "Unavailable after management stops",
            "Applied nourishment",
        ],
        loc="upper left",
        bbox_to_anchor=(1.005, 1.0),
    )
    axes[2].legend(
        handles=[
            Patch(facecolor="#225ea8", label="Both OFF"),
            Patch(facecolor="#41ab5d", label="Both ON"),
            Patch(facecolor="#984ea3", label="Beach ON / Roadway OFF"),
            Patch(facecolor="#ff8c00", label="Beach OFF / Roadway ON"),
            Patch(facecolor="#d9d9d9", label="One or both unavailable"),
        ],
        loc="upper left",
        bbox_to_anchor=(1.005, 1.0),
    )
    figure.suptitle(
        "Dune-migration state: historical-nourishment manager comparison\n"
        "Domains 111–119; yellow circles mark applied nourishment",
        fontsize=14,
    )
    figure.tight_layout(rect=(0, 0, 0.82, 0.95))
    figure.savefig(FIGURE_PATH, dpi=180, bbox_inches="tight")
    plt.close(figure)

    assert FIGURE_PATH.is_file() and FIGURE_PATH.stat().st_size > 0
    assert CSV_PATH.is_file() and CSV_PATH.stat().st_size > 0
