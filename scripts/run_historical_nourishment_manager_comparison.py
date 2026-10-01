#!/usr/bin/env python3
"""Run one side of the historical nourishment manager comparison.

Run this script once with ``--manager beach-dune`` and once with
``--manager roadway``. Each invocation is a separate Python process so module
state and the seeded Barrier3D simulations cannot leak between scenarios.

Both runs retain the established Pea Island management configuration on domains
80--110. The comparison changes only domains 111--119: BeachDuneManager manages
them in one run and the nourishment-enabled RoadwayManager manages them in the
other. Native triggered roadway relocation and established dune rebuilding are
retained. No historical/forced roadway relocation requests are made.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
from pathlib import Path
import sys

import numpy as np
import yaml

SOURCE_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = SOURCE_ROOT.parent
BASE_PROJECT_ROOT = WORKSPACE_ROOT / "CASCADE"
BASE_RUNNER_PATH = (
    BASE_PROJECT_ROOT
    / "scripts"
    / "Pea_Island_ms"
    / "baseline"
    / "run_original_roadway_manager_baseline.py"
)
COMPARISON_ROOT = (
    BASE_PROJECT_ROOT
    / "output"
    / "roadway_nourishment"
    / "historical_manager_comparison"
)

START_YEAR = 1992
END_YEAR = 2007
FIRST_REAL_DOMAIN = 80
START_REAL_INDEX = 71
COMPARISON_DOMAINS = list(range(111, 120))
ORIGINAL_BEACH_DUNE_DOMAINS = [80, 81]
ORIGINAL_ROADWAY_DOMAINS = list(range(82, 111))

SCENARIOS = {
    "beach-dune": {
        "directory": "beach_dune_manager",
        "run_name": (
            "PEA_1992_2007_HistoricalNourishment_" "BeachDuneManager_Hs2p0_Berm1p7"
        ),
        "label": "BeachDuneManager historical nourishment",
    },
    "roadway": {
        "directory": "roadway_manager",
        "run_name": (
            "PEA_1992_2007_HistoricalNourishment_" "RoadwayManager_Hs2p0_Berm1p7"
        ),
        "label": "RoadwayManager historical nourishment",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manager", choices=SCENARIOS, required=True)
    parser.add_argument(
        "--finalize-existing",
        action="store_true",
        help="Validate and finalize a completed NPZ without rerunning it.",
    )
    return parser.parse_args()


def load_base_runner():
    specification = importlib.util.spec_from_file_location(
        "pea_island_historical_nourishment_base_configuration",
        BASE_RUNNER_PATH,
    )
    if specification is None or specification.loader is None:
        raise ImportError(f"Cannot load run configuration: {BASE_RUNNER_PATH}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def real_domain_to_index(domain: int) -> int:
    return START_REAL_INDEX + domain - FIRST_REAL_DOMAIN


def expected_events(base_runner) -> list[dict]:
    events = []
    for domain, years in base_runner.BN_YEARS_BY_DOMAIN.items():
        volumes = base_runner.BN_VOLUME_M3_BY_DOMAIN[domain]
        for year, total_volume_m3 in zip(years, volumes):
            if START_YEAR <= year <= END_YEAR and domain in COMPARISON_DOMAINS:
                events.append(
                    {
                        "year": year,
                        "domain": domain,
                        "total_volume_m3": total_volume_m3,
                        "volume_m3_per_m": total_volume_m3
                        / base_runner.DOMAIN_LENGTH_M,
                    }
                )
    return events


def get_saved_cascade(npz_path: Path):
    source_path = str(SOURCE_ROOT)
    if source_path not in sys.path:
        sys.path.insert(0, source_path)
    with np.load(npz_path, allow_pickle=True) as archive:
        return archive["cascade"].item()


def validate_and_write_events(
    manager_choice: str,
    npz_path: Path,
    output_csv: Path,
    events: list[dict],
) -> dict:
    cascade = get_saved_cascade(npz_path)
    manager_collection = (
        cascade.nourishments if manager_choice == "beach-dune" else cascade.roadways
    )
    rows = []
    unapplied_events = []

    for event in events:
        domain = event["domain"]
        saved_index = real_domain_to_index(domain)
        time_index = event["year"] - START_YEAR + 1
        manager = manager_collection[saved_index]
        event_flags = np.asarray(manager._nourishment_TS, dtype=bool)
        event_volumes = np.asarray(manager._nourishment_volume_TS, dtype=float)
        beach_width = np.asarray(manager.beach_width, dtype=float)
        applied = bool(event_flags[time_index])
        applied_volume = float(event_volumes[time_index])
        volume_matches = bool(np.isclose(applied_volume, event["volume_m3_per_m"]))
        management_stopped = False
        if manager_choice == "roadway":
            management_stopped = bool(cascade.road_break[saved_index])
        else:
            management_stopped = bool(cascade.community_break[saved_index])
        if not applied or not volume_matches:
            unapplied_events.append(
                {
                    "domain": domain,
                    "year": event["year"],
                    "requested_volume_m3_per_m": event["volume_m3_per_m"],
                    "recorded_volume_m3_per_m": applied_volume,
                    "management_stopped_during_run": management_stopped,
                }
            )
        rows.append(
            {
                **event,
                "manager": manager_choice,
                "saved_cascade_index": saved_index,
                "saved_time_index": time_index,
                "nourishment_recorded": applied,
                "recorded_volume_m3_per_m": applied_volume,
                "volume_matches_request": volume_matches,
                "management_stopped_during_run": management_stopped,
                "beach_width_after_event_m": float(beach_width[time_index]),
                "shoreline_after_event_dam": float(
                    cascade.barrier3d[saved_index].x_s_TS[time_index]
                ),
            }
        )

    relocation_events = []
    roadway_mask = np.asarray(cascade.roadway_management_module, dtype=bool)
    for domain in range(80, 120):
        saved_index = real_domain_to_index(domain)
        if roadway_mask[saved_index]:
            manager = cascade.roadways[saved_index]
            for series_name in (
                "triggered_relocation_TS",
                "historical_relocation_requested_TS",
                "forced_relocation_TS",
                "relocation_incomplete_TS",
            ):
                event_indices = np.flatnonzero(
                    np.asarray(getattr(manager, series_name), dtype=bool)
                )
                relocation_events.extend(
                    {
                        "domain": domain,
                        "series": series_name,
                        "time_index": int(index),
                    }
                    for index in event_indices
                )

    with output_csv.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    return {
        "expected_event_count": len(events),
        "applied_event_count": len(events) - len(unapplied_events),
        "unapplied_event_count": len(unapplied_events),
        "unapplied_events": unapplied_events,
        "comparison_domains": COMPARISON_DOMAINS,
        "relocation_event_count": len(relocation_events),
        "relocation_events": relocation_events,
        "event_validation_csv": str(output_csv),
    }


def update_manifest(
    manager_choice: str,
    manifest_path: Path,
    named_manifest_path: Path,
    validation: dict,
) -> None:
    with manifest_path.open() as stream:
        manifest = yaml.safe_load(stream)

    scenario = SCENARIOS[manager_choice]
    identity = manifest["baseline_identity"]
    identity["source_role"] = "cumulative_roadway_functionalities"
    identity["comparison_scenario"] = scenario["label"]
    identity["roadway_manager_relocation_added_used"] = manager_choice == "roadway"
    identity["roadway_nourishment_added_used"] = manager_choice == "roadway"
    identity["cascade_source_sha256"] = sha256(SOURCE_ROOT / "cascade" / "cascade.py")

    simulation = manifest["simulation"]
    simulation["roadway_managed_domains"] = (
        ORIGINAL_ROADWAY_DOMAINS + COMPARISON_DOMAINS
        if manager_choice == "roadway"
        else ORIGINAL_ROADWAY_DOMAINS
    )
    simulation["beach_dune_managed_domains"] = (
        ORIGINAL_BEACH_DUNE_DOMAINS + COMPARISON_DOMAINS
        if manager_choice == "beach-dune"
        else ORIGINAL_BEACH_DUNE_DOMAINS
    )
    simulation["unmanaged_real_domains"] = []
    unapplied_keys = {
        (record["domain"], record["year"]) for record in validation["unapplied_events"]
    }
    historical_requests = []
    for record in simulation["historical_nourishment_requests"]:
        key = (record["domain"], record["year"])
        applied = key not in unapplied_keys
        historical_requests.append(
            {
                **record,
                "applied_by_selected_manager": applied,
                "manager": manager_choice,
                "application_status": (
                    "applied"
                    if applied
                    else "not applied; management stopped during run"
                ),
            }
        )
    simulation["historical_nourishment_requests"] = historical_requests
    simulation["automatic_beach_width_nourishment"] = False
    simulation["nourishment_interval"] = None
    simulation["established_dune_rebuild_behavior_active"] = True
    simulation["historical_or_forced_relocation_requests"] = False
    simulation["native_road_relocation_logic_changed"] = False
    simulation["validation"] = validation

    for record in manifest["initial_model_state_real_domains"]:
        domain = record["domain"]
        if domain in ORIGINAL_BEACH_DUNE_DOMAINS:
            record["manager"] = "BeachDuneManager"
        elif domain in ORIGINAL_ROADWAY_DOMAINS:
            record["manager"] = "RoadwayManager"
        elif manager_choice == "roadway":
            record["manager"] = "RoadwayManager"
        else:
            record["manager"] = "BeachDuneManager"

    with named_manifest_path.open("w") as stream:
        yaml.safe_dump(manifest, stream, sort_keys=False)
    if manifest_path != named_manifest_path:
        manifest_path.unlink()


def main() -> None:
    args = parse_args()
    scenario = SCENARIOS[args.manager]
    output_root = COMPARISON_ROOT / scenario["directory"]
    run_dir = output_root / scenario["run_name"]
    if run_dir.exists() and not args.finalize_existing:
        raise FileExistsError(
            f"Fresh-run protection refused to reuse existing directory: {run_dir}"
        )
    if not run_dir.exists() and args.finalize_existing:
        raise FileNotFoundError(
            f"Cannot finalize because the run directory does not exist: {run_dir}"
        )

    base_runner = load_base_runner()
    manager_path = SOURCE_ROOT / "cascade" / "roadway_manager.py"
    base_runner.ORIGINAL_SOURCE_ROOT = SOURCE_ROOT.resolve()
    base_runner.RUN_NAME = scenario["run_name"]
    base_runner.BASELINE_ROOT = output_root.resolve()
    base_runner.RUN_DIR = run_dir.resolve()
    base_runner.EXPECTED_ROADWAY_MANAGER_SHA256 = base_runner.sha256(manager_path)

    # Preserve management on domains 80--110. Change only which manager applies
    # historical nourishment on comparison domains 111--119.
    if args.manager == "beach-dune":
        base_runner.ROADWAY_MANAGED_DOMAINS = ORIGINAL_ROADWAY_DOMAINS.copy()
        base_runner.BEACH_DUNE_MANAGED_DOMAINS = (
            ORIGINAL_BEACH_DUNE_DOMAINS + COMPARISON_DOMAINS
        )
    else:
        base_runner.ROADWAY_MANAGED_DOMAINS = (
            ORIGINAL_ROADWAY_DOMAINS + COMPARISON_DOMAINS
        )
        base_runner.BEACH_DUNE_MANAGED_DOMAINS = ORIGINAL_BEACH_DUNE_DOMAINS.copy()

    events = expected_events(base_runner)
    if len(events) != 30:
        raise RuntimeError(f"Expected 30 historical events, found {len(events)}")

    if not args.finalize_existing:
        print(f"Comparison scenario: {scenario['label']}")
        print(f"CASCADE source: {SOURCE_ROOT}")
        print(f"Roadway-managed domains: {base_runner.ROADWAY_MANAGED_DOMAINS}")
        print(f"BeachDune-managed domains: {base_runner.BEACH_DUNE_MANAGED_DOMAINS}")
        print("Historical/forced relocation requests: none")
        print("Native triggered relocation: retained")
        print("Established dune-rebuild behavior: retained")
        print(f"Separate output: {run_dir}")
        base_runner.main()
    else:
        print(f"Finalizing existing completed simulation: {run_dir}")

    npz_path = run_dir / f"{scenario['run_name']}.npz"
    validation_csv = run_dir / "historical_nourishment_event_validation.csv"
    validation = validate_and_write_events(
        args.manager,
        npz_path,
        validation_csv,
        events,
    )
    manifest_path = run_dir / "baseline_manifest.yaml"
    named_manifest_path = run_dir / (
        f"historical_nourishment_{scenario['directory']}_manifest.yaml"
    )
    if args.finalize_existing and not manifest_path.is_file():
        manifest_path = named_manifest_path
    update_manifest(
        args.manager,
        manifest_path,
        named_manifest_path,
        validation,
    )
    print(f"Historical nourishment requests: {len(events)}")
    print(f"Applied events: {validation['applied_event_count']}")
    print(f"Unapplied events: {validation['unapplied_event_count']}")
    print(
        "Recorded roadway relocation diagnostic events: "
        f"{validation['relocation_event_count']}"
    )
    print(f"Saved manifest: {named_manifest_path}")
    print(f"Saved event validation: {validation_csv}")


if __name__ == "__main__":
    main()
