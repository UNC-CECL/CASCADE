#!/usr/bin/env python3
"""Plot both historical nourishment simulations and compare overwash."""

from __future__ import annotations

import csv
import importlib.util
from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import yaml

SOURCE_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = SOURCE_ROOT.parent
BASE_PROJECT_ROOT = WORKSPACE_ROOT / "CASCADE"
RESULT_SCRIPT_ROOT = BASE_PROJECT_ROOT / "scripts" / "Pea_Island_ms" / "Result"
SHORELINE_SCRIPT = RESULT_SCRIPT_ROOT / (
    "plot_cascade_shoreline_positions_all_timesteps.py"
)
DUNE_SCRIPT = RESULT_SCRIPT_ROOT / "plot_dunes_one_npz_over_time.py"
OVERWASH_REMOVAL_SCRIPT = RESULT_SCRIPT_ROOT / ("plot_overwash_removal_old_vs_new.py")
COMPARISON_ROOT = (
    BASE_PROJECT_ROOT
    / "output"
    / "roadway_nourishment"
    / "historical_manager_comparison"
)

RUNS = {
    "beach_dune_manager": {
        "name": (
            "PEA_1992_2007_HistoricalNourishment_" "BeachDuneManager_Hs2p0_Berm1p7"
        ),
        "label": "BeachDuneManager on historical nourishment domains",
    },
    "roadway_manager": {
        "name": ("PEA_1992_2007_HistoricalNourishment_" "RoadwayManager_Hs2p0_Berm1p7"),
        "label": "RoadwayManager on historical nourishment domains",
    },
}

FIRST_REAL_INDEX = 71
FIRST_DOMAIN = 80
LAST_DOMAIN = 119
START_YEAR = 1992
END_YEAR = 2007
COMPARISON_DOMAINS = set(range(111, 120))


def load_script(path: Path, module_name: str):
    specification = importlib.util.spec_from_file_location(module_name, path)
    if specification is None or specification.loader is None:
        raise ImportError(f"Cannot load plotting script: {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def run_paths(key: str) -> tuple[Path, Path]:
    run = RUNS[key]
    run_dir = COMPARISON_ROOT / key / run["name"]
    npz_path = run_dir / f"{run['name']}.npz"
    if not npz_path.is_file():
        raise FileNotFoundError(npz_path)
    return run_dir, npz_path


def plot_shorelines() -> None:
    for key, run in RUNS.items():
        run_dir, npz_path = run_paths(key)
        module = load_script(
            SHORELINE_SCRIPT,
            f"plot_cascade_shoreline_{key}",
        )
        module.INPUT_NPZ = npz_path
        module.OUTPUT_DIR = run_dir / "plots" / "shoreline_all_timesteps"
        module.RUN_LABEL = run["label"]
        module.HIGHLIGHT_DOMAINS = tuple(sorted(COMPARISON_DOMAINS))
        module.SHOW_FIGURE = False
        module.main()


def plot_dunes() -> None:
    for key in RUNS:
        run_dir, npz_path = run_paths(key)
        output_dir = run_dir / "plots" / "dune_evolution"
        module = load_script(DUNE_SCRIPT, f"plot_dunes_{key}")
        module.CASCADE_PROJECT_ROOT = SOURCE_ROOT
        module.SHOW_PLOTS = False
        original_argv = sys.argv
        try:
            sys.argv = [
                str(DUNE_SCRIPT),
                "--npz",
                str(npz_path),
                "--output-dir",
                str(output_dir),
                "--domains",
                *[str(domain) for domain in range(FIRST_DOMAIN, LAST_DOMAIN + 1)],
                "--diagnostic-domains",
                *[str(domain) for domain in sorted(COMPARISON_DOMAINS)],
                "--start-year",
                str(START_YEAR),
                "--end-year",
                str(END_YEAR),
                "--gif-fps",
                "0.8",
            ]
            module.main()
        finally:
            sys.argv = original_argv


def plot_overwash_removal() -> None:
    _, beach_npz = run_paths("beach_dune_manager")
    _, roadway_npz = run_paths("roadway_manager")
    output_dir = COMPARISON_ROOT / "comparison_plots" / "overwash_removal"
    module = load_script(
        OVERWASH_REMOVAL_SCRIPT,
        "plot_historical_nourishment_overwash_removal",
    )
    original_argv = sys.argv
    try:
        sys.argv = [
            str(OVERWASH_REMOVAL_SCRIPT),
            "--old-npz",
            str(beach_npz),
            "--new-npz",
            str(roadway_npz),
            "--old-label",
            RUNS["beach_dune_manager"]["label"],
            "--new-label",
            RUNS["roadway_manager"]["label"],
            "--output-dir",
            str(output_dir),
            "--start-year",
            str(START_YEAR),
            "--end-year",
            str(END_YEAR),
        ]
        module.main()
    finally:
        sys.argv = original_argv


def load_cascade(npz_path: Path):
    with np.load(npz_path, allow_pickle=True) as archive:
        return archive["cascade"].item()


def extract_overwash_flux(cascade) -> np.ndarray:
    domains = np.arange(FIRST_DOMAIN, LAST_DOMAIN + 1)
    years = np.arange(START_YEAR, END_YEAR + 1)
    values = np.full((len(domains), len(years)), np.nan, dtype=float)
    for row, domain in enumerate(domains):
        saved_index = FIRST_REAL_INDEX + domain - FIRST_DOMAIN
        series = np.asarray(cascade.barrier3d[saved_index].QowTS, dtype=float)
        count = min(len(years), max(0, len(series) - 1))
        values[row, :count] = series[1 : count + 1]
    return values


def plot_overwash_flux() -> None:
    output_dir = COMPARISON_ROOT / "comparison_plots" / "overwash_flux"
    output_dir.mkdir(parents=True, exist_ok=True)
    _, beach_npz = run_paths("beach_dune_manager")
    _, roadway_npz = run_paths("roadway_manager")
    beach = extract_overwash_flux(load_cascade(beach_npz))
    roadway = extract_overwash_flux(load_cascade(roadway_npz))
    difference = roadway - beach
    domains = np.arange(FIRST_DOMAIN, LAST_DOMAIN + 1)
    years = np.arange(START_YEAR, END_YEAR + 1)

    finite_values = np.concatenate(
        (beach[np.isfinite(beach)], roadway[np.isfinite(roadway)])
    )
    raw_min = float(np.min(finite_values))
    raw_max = float(np.max(finite_values))
    difference_limit = float(np.nanmax(np.abs(difference)))
    if difference_limit == 0.0:
        difference_limit = 1.0

    fig, axes = plt.subplots(3, 1, figsize=(17, 14), layout="constrained")
    extent = (
        years[0] - 0.5,
        years[-1] + 0.5,
        domains[0] - 0.5,
        domains[-1] + 0.5,
    )
    for ax, values, title in zip(
        axes[:2],
        (beach, roadway),
        (
            RUNS["beach_dune_manager"]["label"],
            RUNS["roadway_manager"]["label"],
        ),
    ):
        image = ax.imshow(
            values,
            origin="lower",
            aspect="auto",
            interpolation="nearest",
            extent=extent,
            cmap="viridis",
            vmin=raw_min,
            vmax=raw_max,
        )
        ax.set_title(title)
        fig.colorbar(image, ax=ax, label="Net Barrier3D overwash flux, Qow (m³/m)")

    difference_image = axes[2].imshow(
        difference,
        origin="lower",
        aspect="auto",
        interpolation="nearest",
        extent=extent,
        cmap="RdBu_r",
        norm=TwoSlopeNorm(
            vmin=-difference_limit,
            vcenter=0.0,
            vmax=difference_limit,
        ),
    )
    axes[2].set_title("Difference: RoadwayManager − BeachDuneManager")
    fig.colorbar(
        difference_image,
        ax=axes[2],
        label="Qow difference (m³/m)",
    )
    for ax in axes:
        ax.set_ylabel("Real Pea Island domain")
        ax.set_yticks(domains[::2])
        ax.set_xticks(years)
        ax.set_xticklabels(years, rotation=45, ha="right")
        for boundary in (110.5,):
            ax.axhline(boundary, color="white", linewidth=1.5, linestyle="--")
    axes[-1].set_xlabel("Annual model update")
    fig.suptitle(
        "Historical nourishment manager comparison: net overwash flux",
        fontsize=17,
        fontweight="bold",
    )
    heatmap_path = output_dir / "overwash_flux_heatmaps_and_difference.png"
    fig.savefig(heatmap_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    target_rows = np.array([domain in COMPARISON_DOMAINS for domain in domains])
    fig, axes = plt.subplots(2, 1, figsize=(15, 10), layout="constrained")
    for values, label, color in (
        (beach, RUNS["beach_dune_manager"]["label"], "#4c72b0"),
        (roadway, RUNS["roadway_manager"]["label"], "#dd8452"),
    ):
        axes[0].plot(
            years,
            np.nansum(values, axis=0),
            marker="o",
            linewidth=2,
            label=label,
            color=color,
        )
        axes[1].plot(
            years,
            np.nansum(values[target_rows], axis=0),
            marker="o",
            linewidth=2,
            label=label,
            color=color,
        )
    axes[0].set_title("All real domains 80–119")
    axes[1].set_title("Historical nourishment domains 111–119")
    for ax in axes:
        ax.set_ylabel("Annual summed Qow (m³/m)")
        ax.set_xlabel("Annual model update")
        ax.set_xticks(years)
        ax.set_xticklabels(years, rotation=45, ha="right")
        ax.grid(True, alpha=0.3)
        ax.legend(frameon=False)
    summary_path = output_dir / "overwash_flux_annual_comparison.png"
    fig.savefig(summary_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    csv_path = output_dir / "overwash_flux_domain_year_comparison.csv"
    with csv_path.open("w", newline="") as stream:
        fieldnames = (
            "year",
            "domain_id",
            "historical_nourishment_domain",
            "beach_dune_manager_Qow_m3_per_m",
            "roadway_manager_Qow_m3_per_m",
            "roadway_minus_beach_dune_Qow_m3_per_m",
        )
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for row, domain in enumerate(domains):
            for column, year in enumerate(years):
                writer.writerow(
                    {
                        "year": int(year),
                        "domain_id": int(domain),
                        "historical_nourishment_domain": bool(
                            domain in COMPARISON_DOMAINS
                        ),
                        "beach_dune_manager_Qow_m3_per_m": float(beach[row, column]),
                        "roadway_manager_Qow_m3_per_m": float(roadway[row, column]),
                        "roadway_minus_beach_dune_Qow_m3_per_m": float(
                            difference[row, column]
                        ),
                    }
                )

    metrics = {
        "quantity": "Barrier3D QowTS (net overwash flux after management)",
        "units": "m3/m",
        "all_domains_80_119": {
            "beach_dune_manager_sum": float(np.nansum(beach)),
            "roadway_manager_sum": float(np.nansum(roadway)),
            "roadway_minus_beach_dune_sum": float(np.nansum(difference)),
        },
        "comparison_domains_111_119": {
            "beach_dune_manager_sum": float(np.nansum(beach[target_rows])),
            "roadway_manager_sum": float(np.nansum(roadway[target_rows])),
            "roadway_minus_beach_dune_sum": float(np.nansum(difference[target_rows])),
        },
    }
    metrics_path = output_dir / "overwash_flux_summary.yaml"
    with metrics_path.open("w") as stream:
        yaml.safe_dump(metrics, stream, sort_keys=False)

    print("Saved net-overwash comparison:")
    for path in (heatmap_path, summary_path, csv_path, metrics_path):
        print(f"  {path}")


def main() -> None:
    source_path = str(SOURCE_ROOT)
    if source_path not in sys.path:
        sys.path.insert(0, source_path)
    plot_shorelines()
    plot_dunes()
    plot_overwash_removal()
    plot_overwash_flux()


if __name__ == "__main__":
    main()
