"""Compare roadway and BeachDuneManager nourishment for historical events.

These tests isolate the nourishment operation. Both managers receive identical
pre-event barrier states and the same historical volume in m^3/m. Roadway
relocation requests, BeachDuneManager overwash filtering, and dune rebuilding
are intentionally not invoked.
"""

from copy import deepcopy
import csv
from pathlib import Path
from types import SimpleNamespace

import matplotlib
import numpy as np
import pytest
from matplotlib.colors import ListedColormap
from matplotlib.patches import FancyBboxPatch

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from cascade.beach_dune_manager import BeachDuneManager
from cascade.roadway_manager import RoadwayManager

DOMAIN_LENGTH_M = 500.0
WORKSPACE_ROOT = Path(__file__).resolve().parents[2]
CONTROLLED_TEST_OUTPUT = (
    WORKSPACE_ROOT
    / "CASCADE"
    / "output"
    / "roadway_nourishment"
    / "controlled_nourishment_test"
)

# Total project volumes and years from the Pea Island historical nourishment
# record. Events after the 1992-2007 hindcast are excluded from this test.
HISTORICAL_NOURISHMENT = {
    111: [(2003, 102954.25)],
    112: [(2003, 102954.25)],
    113: [(2002, 96763.6), (2003, 102954.25)],
    114: [(1992, 431200.0), (2002, 96763.6), (2003, 102954.25)],
    115: [
        (1992, 431200.0),
        (2002, 96763.6),
        (2003, 102954.25),
        (2004, 74972.4),
    ],
    116: [
        (1992, 49146.6667),
        (1993, 173294.0),
        (2001, 102741.25),
        (2002, 96763.6),
        (2003, 102954.25),
        (2004, 74972.4),
    ],
    117: [
        (1992, 49146.6667),
        (1993, 173294.0),
        (1995, 52184.8),
        (2001, 102741.25),
        (2002, 96763.6),
        (2003, 102954.25),
        (2004, 74972.4),
    ],
    118: [
        (1992, 49146.6667),
        (2001, 102741.25),
        (2003, 102954.25),
        (2004, 74972.4),
    ],
    119: [(2001, 102741.25), (2004, 74972.4)],
}

HISTORICAL_EVENTS = [
    pytest.param(
        domain,
        year,
        total_volume_m3,
        id=f"domain-{domain}-year-{year}",
    )
    for domain, events in HISTORICAL_NOURISHMENT.items()
    for year, total_volume_m3 in events
]


def nourishment_test_barrier():
    """Return the minimum identical post-storm state needed by both managers."""

    time_step_count = 4
    alongshore_cells = 5
    interior = np.full((12, alongshore_cells), 0.2)
    dunes = np.full((time_step_count, alongshore_cells, 2), 0.4)
    domain_ts = np.empty(time_step_count, dtype=object)
    for index in range(time_step_count):
        domain_ts[index] = interior.copy()

    return SimpleNamespace(
        time_index=2,
        x_s=100.0,
        x_t=0.0,
        x_s_TS=[100.0, 100.0],
        x_b_TS=[114.0, 114.0],
        s_sf_TS=[0.01, 0.01],
        h_b_TS=[0.2, 0.2],
        InteriorWidth_AvgTS=[12.0, 12.0],
        DShoreface=1.0,
        DuneDomain=dunes,
        InteriorDomain=interior,
        DomainTS=domain_ts,
        growthparam=np.full((1, alongshore_cells), 0.5),
        QowTS=[0.0, 0.0],
        dune_migration_on=True,
        SCRagg=np.zeros(time_step_count),
    )


@pytest.mark.parametrize("domain,year,total_volume_m3", HISTORICAL_EVENTS)
def test_historical_nourishment_event_matches_beach_dune_manager(
    domain,
    year,
    total_volume_m3,
):
    """Both managers must apply each historical nourishment event identically."""

    del domain, year  # Values identify each independently reported pytest case.
    nourishment_volume_m3_per_m = total_volume_m3 / DOMAIN_LENGTH_M

    beach_barrier = nourishment_test_barrier()
    roadway_barrier = deepcopy(beach_barrier)

    beach_manager = BeachDuneManager(
        nourishment_interval=None,
        nourishment_volume=nourishment_volume_m3_per_m,
        initial_beach_width=30.0,
        time_step_count=4,
    )
    # Disable BeachDuneManager processes that are not part of nourishment.
    beach_manager.overwash_removal = False

    roadway_manager = RoadwayManager(
        road_width=10.0,
        road_setback=5.0,
        nourishment_interval=None,
        nourishment_volume=nourishment_volume_m3_per_m,
        initial_beach_width=30.0,
        time_step_count=4,
    )

    beach_request, rebuild_request = beach_manager.update(
        barrier3d=beach_barrier,
        nourish_now=True,
        rebuild_dune_now=False,
        nourishment_interval=None,
    )
    # Call only the roadway nourishment operation; relocation logic is not run.
    roadway_manager.nourish_now(roadway_barrier)

    event_index = beach_barrier.time_index - 1
    assert beach_request == 0
    assert rebuild_request is False
    assert beach_manager._nourishment_TS[event_index]
    assert roadway_manager.nourishment_TS[event_index]
    assert beach_manager.nourishment_volume_TS[event_index] == pytest.approx(
        nourishment_volume_m3_per_m
    )
    assert roadway_manager.nourishment_volume_TS[event_index] == pytest.approx(
        nourishment_volume_m3_per_m
    )
    assert roadway_barrier.x_s == pytest.approx(beach_barrier.x_s)
    assert roadway_barrier.x_s_TS[-1] == pytest.approx(beach_barrier.x_s_TS[-1])
    assert roadway_barrier.s_sf_TS[-1] == pytest.approx(beach_barrier.s_sf_TS[-1])
    assert roadway_manager.beach_width[event_index] == pytest.approx(
        beach_manager.beach_width[event_index]
    )
    assert roadway_barrier.dune_migration_on is False
    assert beach_barrier.dune_migration_on is False
    assert not roadway_manager.historical_relocation_requested_TS.any()
    assert not roadway_manager.forced_relocation_TS.any()
    assert not roadway_manager.triggered_relocation_TS.any()
    assert not roadway_manager.relocation_incomplete_TS.any()


def test_historical_nourishment_schedule_covers_1992_through_2007_only():
    """Guard the event count and hindcast-year boundary used by the comparison."""

    events = [
        event
        for domain_events in HISTORICAL_NOURISHMENT.values()
        for event in domain_events
    ]

    assert len(events) == 30
    assert all(1992 <= year <= 2007 for year, _ in events)


def run_controlled_nourishment_comparison(total_volume_m3):
    """Run the exact controlled A/B operation and return its recorded outputs."""

    nourishment_volume_m3_per_m = total_volume_m3 / DOMAIN_LENGTH_M
    beach_barrier = nourishment_test_barrier()
    roadway_barrier = deepcopy(beach_barrier)

    beach_manager = BeachDuneManager(
        nourishment_interval=None,
        nourishment_volume=nourishment_volume_m3_per_m,
        initial_beach_width=30.0,
        time_step_count=4,
    )
    beach_manager.overwash_removal = False
    roadway_manager = RoadwayManager(
        road_width=10.0,
        road_setback=5.0,
        nourishment_interval=None,
        nourishment_volume=nourishment_volume_m3_per_m,
        initial_beach_width=30.0,
        time_step_count=4,
    )

    beach_manager.update(
        barrier3d=beach_barrier,
        nourish_now=True,
        rebuild_dune_now=False,
        nourishment_interval=None,
    )
    roadway_manager.nourish_now(roadway_barrier)
    event_index = beach_barrier.time_index - 1

    return {
        "input_volume_m3_per_m": nourishment_volume_m3_per_m,
        "beach_x_s_dam": float(beach_barrier.x_s),
        "roadway_x_s_dam": float(roadway_barrier.x_s),
        "beach_s_sf": float(beach_barrier.s_sf_TS[-1]),
        "roadway_s_sf": float(roadway_barrier.s_sf_TS[-1]),
        "beach_width_m": float(beach_manager.beach_width[event_index]),
        "roadway_width_m": float(roadway_manager.beach_width[event_index]),
        "beach_applied_volume_m3_per_m": float(
            beach_manager.nourishment_volume_TS[event_index]
        ),
        "roadway_applied_volume_m3_per_m": float(
            roadway_manager.nourishment_volume_TS[event_index]
        ),
        "beach_event_flag": bool(beach_manager._nourishment_TS[event_index]),
        "roadway_event_flag": bool(roadway_manager.nourishment_TS[event_index]),
        "beach_dune_migration_on": bool(beach_barrier.dune_migration_on),
        "roadway_dune_migration_on": bool(roadway_barrier.dune_migration_on),
        "roadway_relocation_diagnostics_false": bool(
            not roadway_manager.historical_relocation_requested_TS.any()
            and not roadway_manager.forced_relocation_TS.any()
            and not roadway_manager.triggered_relocation_TS.any()
            and not roadway_manager.relocation_incomplete_TS.any()
        ),
    }


def draw_workflow_box(axis, x, y, width, height, text, color):
    box = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.015",
        linewidth=1.5,
        edgecolor="#333333",
        facecolor=color,
        transform=axis.transAxes,
    )
    axis.add_patch(box)
    axis.text(
        x + width / 2,
        y + height / 2,
        text,
        transform=axis.transAxes,
        ha="center",
        va="center",
        fontsize=10,
    )


def create_controlled_test_design_figure(example):
    barrier = nourishment_test_barrier()
    interior_m = barrier.InteriorDomain * 10.0
    dunes_m = barrier.DuneDomain[barrier.time_index - 1] * 10.0

    figure = plt.figure(figsize=(17, 10))
    grid = figure.add_gridspec(2, 3, height_ratios=[1.0, 1.1], hspace=0.42)

    interior_axis = figure.add_subplot(grid[0, 0])
    interior_image = interior_axis.imshow(
        interior_m.T,
        aspect="auto",
        cmap="terrain",
        vmin=0,
        vmax=5,
        interpolation="none",
    )
    interior_axis.set_title("Synthetic interior grid: 12 × 5 cells")
    interior_axis.set_xlabel("Cross-shore cell")
    interior_axis.set_ylabel("Alongshore cell")
    figure.colorbar(
        interior_image,
        ax=interior_axis,
        label="Elevation (m MHW)",
        shrink=0.85,
    )
    interior_axis.text(
        0.5,
        -0.23,
        "Every cell = 0.2 dam = 2.0 m MHW",
        transform=interior_axis.transAxes,
        ha="center",
        fontsize=10,
    )

    dune_axis = figure.add_subplot(grid[0, 1])
    dune_image = dune_axis.imshow(
        dunes_m,
        aspect="auto",
        cmap="YlOrBr",
        vmin=0,
        vmax=5,
        interpolation="none",
    )
    dune_axis.set_title("Synthetic dune grid: 5 × 2 cells")
    dune_axis.set_xlabel("Dune row")
    dune_axis.set_ylabel("Alongshore cell")
    dune_axis.set_xticks([0, 1], ["Front", "Back"])
    figure.colorbar(
        dune_image,
        ax=dune_axis,
        label="Height above berm (m)",
        shrink=0.85,
    )
    dune_axis.text(
        0.5,
        -0.23,
        "Every cell = 0.4 dam = 4.0 m above berm",
        transform=dune_axis.transAxes,
        ha="center",
        fontsize=10,
    )

    parameter_axis = figure.add_subplot(grid[0, 2])
    parameter_axis.axis("off")
    parameter_axis.set_title("Identical starting values")
    parameter_text = (
        "Synthetic state — not a real Pea Island domain\n\n"
        "Shoreline, xₛ: 100 dam\n"
        "Shoreface toe, xₜ: 0 dam\n"
        "Shoreface slope: 0.010\n"
        "Shoreface depth: 1 dam\n"
        "Average barrier height: 0.2 dam\n"
        "Beach width: 30 m\n"
        "Dune migration: ON\n\n"
        "Example historical volume:\n"
        "Domain 114, 1992\n"
        "431,200 m³ ÷ 500 m = 862.4 m³/m"
    )
    parameter_axis.text(
        0.03,
        0.93,
        parameter_text,
        va="top",
        fontsize=11,
        linespacing=1.35,
        bbox={"boxstyle": "round,pad=0.6", "facecolor": "#f7f7f7"},
    )

    workflow_axis = figure.add_subplot(grid[1, :])
    workflow_axis.axis("off")
    draw_workflow_box(
        workflow_axis,
        0.01,
        0.35,
        0.17,
        0.3,
        "ONE SYNTHETIC\nBARRIER STATE",
        "#deebf7",
    )
    draw_workflow_box(
        workflow_axis,
        0.25,
        0.35,
        0.14,
        0.3,
        "DEEP COPY\nidentical inputs",
        "#f7f7f7",
    )
    draw_workflow_box(
        workflow_axis,
        0.47,
        0.60,
        0.21,
        0.25,
        "Original BeachDuneManager\nnourish_now=True\noverwash OFF; rebuild OFF",
        "#c7e9c0",
    )
    draw_workflow_box(
        workflow_axis,
        0.47,
        0.15,
        0.21,
        0.25,
        "Modified RoadwayManager\nnourish_now() only\nroad processes bypassed",
        "#fdd0a2",
    )
    output_text = (
        "MATCHED OUTPUTS\n"
        f"xₛ = {example['beach_x_s_dam']:.5f} dam\n"
        f"s_sf = {example['beach_s_sf']:.9f}\n"
        f"beach width = {example['beach_width_m']:.2f} m\n"
        f"applied = {example['input_volume_m3_per_m']:.1f} m³/m\n"
        "event flag = True\n"
        "dune migration = OFF\n"
        "relocation flags = False"
    )
    draw_workflow_box(
        workflow_axis,
        0.77,
        0.27,
        0.21,
        0.46,
        output_text,
        "#d9f0d3",
    )
    arrow = {"arrowstyle": "->", "linewidth": 2, "color": "#444444"}
    workflow_axis.annotate(
        "", xy=(0.25, 0.5), xytext=(0.18, 0.5), xycoords="axes fraction", arrowprops=arrow
    )
    workflow_axis.annotate(
        "",
        xy=(0.47, 0.725),
        xytext=(0.39, 0.53),
        xycoords="axes fraction",
        arrowprops=arrow,
    )
    workflow_axis.annotate(
        "",
        xy=(0.47, 0.275),
        xytext=(0.39, 0.47),
        xycoords="axes fraction",
        arrowprops=arrow,
    )
    workflow_axis.annotate(
        "",
        xy=(0.77, 0.57),
        xytext=(0.68, 0.725),
        xycoords="axes fraction",
        arrowprops=arrow,
    )
    workflow_axis.annotate(
        "",
        xy=(0.77, 0.43),
        xytext=(0.68, 0.275),
        xycoords="axes fraction",
        arrowprops=arrow,
    )
    workflow_axis.text(
        0.575,
        0.48,
        "same 862.4 m³/m",
        transform=workflow_axis.transAxes,
        ha="center",
        fontsize=10,
        fontweight="bold",
    )
    figure.suptitle(
        "Controlled A/B test of the nourishment operation\n"
        "Both managers receive the same synthetic barrier and the same sand volume",
        fontsize=16,
        fontweight="bold",
    )
    return figure


def create_all_events_result_figure(records):
    event_labels = [f"D{record['domain']}\n{record['year']}" for record in records]
    x = np.arange(len(records))
    volumes = np.asarray([record["input_volume_m3_per_m"] for record in records])

    figure = plt.figure(figsize=(18, 12))
    grid = figure.add_gridspec(3, 2, height_ratios=[1, 1, 1.15], hspace=0.38)
    volume_axis = figure.add_subplot(grid[0, 0])
    shoreline_axis = figure.add_subplot(grid[0, 1])
    width_axis = figure.add_subplot(grid[1, 0])
    slope_axis = figure.add_subplot(grid[1, 1])
    pass_axis = figure.add_subplot(grid[2, :])

    volume_axis.bar(x, volumes, color="#6baed6")
    volume_axis.set_title("Historical volume supplied to both managers")
    volume_axis.set_ylabel("Nourishment (m³/m)")
    volume_axis.set_xticks([])

    comparisons = [
        (
            shoreline_axis,
            "beach_x_s_dam",
            "roadway_x_s_dam",
            "Shoreline after nourishment",
            "xₛ (dam)",
        ),
        (
            width_axis,
            "beach_width_m",
            "roadway_width_m",
            "Beach width after nourishment",
            "Beach width (m)",
        ),
        (
            slope_axis,
            "beach_s_sf",
            "roadway_s_sf",
            "Shoreface slope after nourishment",
            "Slope (dimensionless)",
        ),
    ]
    for axis, beach_key, roadway_key, title, ylabel in comparisons:
        beach_values = [record[beach_key] for record in records]
        roadway_values = [record[roadway_key] for record in records]
        axis.plot(x, beach_values, "o-", color="#238b45", label="BeachDuneManager")
        axis.plot(
            x,
            roadway_values,
            "x",
            color="#d95f0e",
            markersize=8,
            markeredgewidth=1.5,
            label="RoadwayManager",
        )
        axis.set_title(title)
        axis.set_ylabel(ylabel)
        axis.set_xticks([])
        axis.legend()

    pass_names = [
        "Shoreline position",
        "Shoreface slope",
        "Beach width",
        "Applied volume",
        "Event flag",
        "Dune migration OFF",
        "Relocation diagnostics false",
    ]
    pass_matrix = np.ones((len(pass_names), len(records)))
    pass_axis.imshow(
        pass_matrix,
        aspect="auto",
        interpolation="none",
        cmap=ListedColormap(["#31a354"]),
        vmin=0,
        vmax=1,
    )
    pass_axis.set_title("Assertion result for every event: all green cells PASSED")
    pass_axis.set_yticks(np.arange(len(pass_names)), pass_names)
    pass_axis.set_xticks(x, event_labels, rotation=90)
    pass_axis.set_xlabel("Historical event identifier (domain and year)")
    for row in range(len(pass_names)):
        for column in range(len(records)):
            pass_axis.text(
                column,
                row,
                "✓",
                ha="center",
                va="center",
                color="white",
                fontsize=8,
                fontweight="bold",
            )

    figure.suptitle(
        "Controlled nourishment-function comparison: all 30 historical volumes\n"
        "Overlapping markers mean both implementations returned the same value",
        fontsize=16,
        fontweight="bold",
    )
    figure.subplots_adjust(left=0.09, right=0.97, bottom=0.12, top=0.91)
    return figure


def test_visualize_controlled_nourishment_comparison():
    """Create advisor-ready documentation from the exact synthetic A/B test."""

    records = []
    for domain, events in HISTORICAL_NOURISHMENT.items():
        for year, total_volume_m3 in events:
            record = run_controlled_nourishment_comparison(total_volume_m3)
            record.update(
                {
                    "domain": domain,
                    "year": year,
                    "total_volume_m3": total_volume_m3,
                }
            )

            assert record["roadway_x_s_dam"] == pytest.approx(
                record["beach_x_s_dam"]
            )
            assert record["roadway_s_sf"] == pytest.approx(record["beach_s_sf"])
            assert record["roadway_width_m"] == pytest.approx(record["beach_width_m"])
            assert record["roadway_applied_volume_m3_per_m"] == pytest.approx(
                record["beach_applied_volume_m3_per_m"]
            )
            assert record["beach_event_flag"] and record["roadway_event_flag"]
            assert not record["beach_dune_migration_on"]
            assert not record["roadway_dune_migration_on"]
            assert record["roadway_relocation_diagnostics_false"]
            records.append(record)

    assert len(records) == 30
    CONTROLLED_TEST_OUTPUT.mkdir(parents=True, exist_ok=True)
    design_path = CONTROLLED_TEST_OUTPUT / "controlled_test_design_and_example.png"
    results_path = CONTROLLED_TEST_OUTPUT / "all_30_historical_volumes_results.png"
    csv_path = CONTROLLED_TEST_OUTPUT / "all_30_historical_volumes_results.csv"

    design_figure = create_controlled_test_design_figure(records[4])
    design_figure.savefig(design_path, dpi=190, bbox_inches="tight")
    plt.close(design_figure)

    results_figure = create_all_events_result_figure(records)
    results_figure.savefig(results_path, dpi=190, bbox_inches="tight")
    plt.close(results_figure)

    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0].keys()))
        writer.writeheader()
        writer.writerows(records)

    assert design_path.is_file() and design_path.stat().st_size > 0
    assert results_path.is_file() and results_path.stat().st_size > 0
    assert csv_path.is_file() and csv_path.stat().st_size > 0
