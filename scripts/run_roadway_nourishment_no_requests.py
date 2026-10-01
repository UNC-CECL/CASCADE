#!/usr/bin/env python3
"""Run the current RoadwayManager source with no nourishment requests.

This is the regression control for the relocation plus roadway-nourishment source.
It uses the established original-roadway initial conditions, explicitly supplies no
historical roadway nourishment requests, and writes to a new, separate directory.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

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
RUN_NAME = "PEA_1992_2007_RoadwayIndependentNourishment_NoRequests_Hs2p0_Berm1p7"
OUTPUT_ROOT = BASE_PROJECT_ROOT / "output" / "roadway_nourishment" / "no_requests"
RUN_DIR = OUTPUT_ROOT / RUN_NAME


def load_base_runner():
    specification = importlib.util.spec_from_file_location(
        "pea_island_original_baseline_configuration",
        BASE_RUNNER_PATH,
    )
    if specification is None or specification.loader is None:
        raise ImportError(f"Cannot load run configuration: {BASE_RUNNER_PATH}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def main():
    manager_path = SOURCE_ROOT / "cascade" / "roadway_manager.py"
    if not manager_path.is_file():
        raise FileNotFoundError(manager_path)

    runner = load_base_runner()
    runner.ORIGINAL_SOURCE_ROOT = SOURCE_ROOT.resolve()
    runner.RUN_NAME = RUN_NAME
    runner.BASELINE_ROOT = OUTPUT_ROOT.resolve()
    runner.RUN_DIR = RUN_DIR.resolve()
    runner.EXPECTED_ROADWAY_MANAGER_SHA256 = runner.sha256(manager_path)

    # This control contains no manual nourishment requests. The original baseline
    # assigned requests to roadway domains, where the original manager ignored
    # them; the enhanced manager would correctly act on them.
    runner.BN_YEARS_BY_DOMAIN = {}
    runner.BN_VOLUME_M3_BY_DOMAIN = {}

    print("Current RoadwayManager control: no nourishment requests")
    print(f"CASCADE source: {runner.ORIGINAL_SOURCE_ROOT}")
    print(f"Separate output: {runner.RUN_DIR}")
    runner.main()

    manifest_path = RUN_DIR / "baseline_manifest.yaml"
    with manifest_path.open() as stream:
        manifest = yaml.safe_load(stream)
    identity = manifest["baseline_identity"]
    identity["source_role"] = "relocation_plus_roadway_nourishment"
    identity["roadway_manager_relocation_added_used"] = True
    identity["roadway_nourishment_added_used"] = True
    simulation = manifest["simulation"]
    simulation["historical_nourishment_requests"] = []
    simulation["automatic_beach_width_nourishment"] = False
    simulation["roadway_nourishment_interval"] = None
    with manifest_path.open("w") as stream:
        yaml.safe_dump(manifest, stream, sort_keys=False)
    named_manifest_path = RUN_DIR / "roadway_nourishment_no_requests_manifest.yaml"
    manifest_path.rename(named_manifest_path)
    print(f"Saved current-source manifest: {named_manifest_path}")


if __name__ == "__main__":
    main()
