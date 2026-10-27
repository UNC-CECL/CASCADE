#!/usr/bin/env python3
"""
A Hatteras hindcast run's checks and written outputs, once the model is built and once it has run.

    from HAT_hindcast_outputs import pre_run_diagnostics, verify_run, write_run_outputs

Each function is one section of HAT_hindcast.py (11.3, 12.2-12.3, 12.4) and works on the run's `r`
(HAT_hindcast_setup.prepare_run, plus the runner's RUN_NAME, RUN_DIR, cascade and timing). Moved out
of HAT_hindcast.py on 2026-10-23, unchanged. Details: scripts/real-world-cases/Hatteras/hatteras_ms/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-23
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd

import cascade.beach_dune_manager as _bdm_module  # dune-cap provenance

from cascade_pipeline import nourishment, roadway
from cascade_pipeline import reports
from cascade_pipeline.hindcast import (
    DAM_TO_M,
    USE_SANDBOX_CASCADE,
    brie_r_ipl,
    measure_groin_extent,
)
from cascade_pipeline.run_layout import write_path
from cascade_pipeline.run_registry import (
    INDEX_SECTION,
    git_provenance,
    skill_vs_target,
    values_digest,
    timestamp,
    write_run_metadata,
)
from cascade_pipeline.shoreline import (build_shoreline_matrix,
                                        compute_change_rate, compute_lrr)

from site_layer.hatteras_site_config import (
    HATTERAS_BE_EDGE_DOMAINS,
    HATTERAS_DOMAINS,
    HATTERAS_FIRST_ROAD_DOMAIN,
    HATTERAS_LAST_ROAD_DOMAIN,
    HATTERAS_GEOMETRY,
    SCORE_INTERIOR_GIS,
)
from run_figures import OVERWASH_NPY, RELOCATIONS_NPY

# The git repository, for the provenance stamp; the scripts are tracked there
REPO_DIR = next(_p for _p in Path(__file__).resolve().parents if (_p / "pyproject.toml").exists())


# Section 11.3: read-outs of the built model, and the groin prediction written down before the run
def pre_run_diagnostics(r):
    # 11.3 Pre-run diagnostics, and the groin prediction

    # Read off the built model; the fillet prediction is written down before the run

    # r_ipl at the groin cell's starting angle (shore-normal is negative under option A)
    if r.GROIN_ENABLED:
        _up = HATTERAS_DOMAINS.gis_to_pad(r.GROIN_UPDRIFT_GIS)
        r.R_IPL_THETA_DEG = float(np.degrees(np.arctan2(
            r.island_offset[_up + 1] - r.island_offset[_up],
            HATTERAS_DOMAINS.domain_spacing_m)))
    else:
        r.R_IPL_THETA_DEG = 0.0
    r.R_IPL = brie_r_ipl(r.cascade, theta_deg=r.R_IPL_THETA_DEG)
    _brie = r.cascade._brie_coupler._brie
    _d_sf_m = float(_brie.d_sf)
    _h_b_m = float(r.cascade.barrier3d[0].h_b_TS[0]) * DAM_TO_M
    _profile_height_m = _d_sf_m + _h_b_m
    _berm_floor_m = float(r.cascade.barrier3d[0].BermEl) * DAM_TO_M

    reports.pre_run_report(
        run_name=r.RUN_NAME, run_dir=r.RUN_DIR, run_years=r.RUN_YEARS,
        start_year=r.START_YEAR, end_year=r.END_YEAR, geometry=HATTERAS_DOMAINS,
        roadway_on=r.ROADWAY_MANAGEMENT_ON, beach_dune_on=r.BEACH_DUNE_MANAGEMENT_ON,
        wave_height=r.Hs, d_sf_m=_d_sf_m, h_b_m=_h_b_m,
        profile_height_m=_profile_height_m, berm_floor_m=_berm_floor_m,
        dune_design_elevation_m=r.DUNE_DESIGN_ELEVATION_M,
        dune_minimum_elevation_m=r.DUNE_MINIMUM_ELEVATION_M,
        groin=r.cascade.groin, groin_enabled=r.GROIN_ENABLED, r_ipl=r.R_IPL)


# Sections 12.2-12.3: the run's own checks, and the groin's extent against its no-groin twin
def verify_run(r):
    # 12.2 Verify

    r.shoreline_m = build_shoreline_matrix(r.cascade)
    r.states, _ = r.shoreline_m.shape

    reports.run_length_report(states=r.states, run_years=r.RUN_YEARS)

    # The denominator behind the old 5% bias: checked, not trusted
    assert r.states - 1 == r.RUN_YEARS, (
        f"span mismatch: {r.states} states implies {r.states - 1} years elapsed, "
        f"but RUN_YEARS is {r.RUN_YEARS}")
    # Two estimators: change_rate (endpoint, conserves net movement) and model_lrr (OLS, matches the target)
    r.change_rate = compute_change_rate(
        r.shoreline_m, span_years=r.RUN_YEARS, flip_sign=r.FLIP_SIGN_MODEL)
    r.model_lrr, r.model_lrr_r2 = compute_lrr(
        r.shoreline_m, span_years=r.RUN_YEARS, flip_sign=r.FLIP_SIGN_MODEL)

    # Nourishment: the model's own record, not the schedule's intent
    r.BN_REPORT = nourishment.verify_nourishment(
        r.cascade, r.BN_SCHEDULE_APPLIED, r.BEACH_DUNE_MANAGEMENT_ON,
        roadway_on=r.ROADWAY_MANAGEMENT_ON)
    reports.nourishment_report(report=r.BN_REPORT, run_years=r.RUN_YEARS)
    assert r.BN_REPORT["ok"], "nourishment did not reach the model as scheduled"

    # The double-management consequence section 6 predicted
    reports.frozen_setbacks_report(
        double_managed=r.DOUBLE_MANAGED_GIS,
        rows=(nourishment.verify_setbacks_frozen(
                  r.cascade, r.DOUBLE_MANAGED_GIS, HATTERAS_DOMAINS)
              if r.DOUBLE_MANAGED_GIS else ()))

    # Which road_offset survived
    r.ROAD_SUMMARY = roadway.summarise_road_management(
        r.cascade, HATTERAS_DOMAINS, HATTERAS_FIRST_ROAD_DOMAIN,
        HATTERAS_LAST_ROAD_DOMAIN)
    r.drowned, r.blocked = reports.roadway_outcome_report(summary=r.ROAD_SUMMARY)

    # 12.3 Groin: the pre-registered extent check

    GROIN_EXTENT_THRESHOLD_FRAC = 0.10   # fraction of peak effect defining "extent"

    r.GROIN_EXTENT = None
    if r.cascade.groin is None:
        print(f"\nGROIN EXTENT          skipped: groin not enabled in this run")
    elif not r.GROIN_BASELINE_NPY:
        print(f"\nGROIN EXTENT          skipped: no paired no-groin baseline")
        print(f"  expected            {r.GROIN_BASELINE_NAME}")
        print(f"  Run once with GROIN_ENABLED = False to create it, then re-run.")
    else:
        _baseline_m = np.load(r.GROIN_BASELINE_NPY)
        r.GROIN_EXTENT = measure_groin_extent(
            r.shoreline_m, _baseline_m, HATTERAS_DOMAINS,
            r.GROIN_UPDRIFT_GIS, r.GROIN_DOWNDRIFT_GIS, GROIN_EXTENT_THRESHOLD_FRAC)
        reports.groin_extent_report(
            extent=r.GROIN_EXTENT, threshold_frac=GROIN_EXTENT_THRESHOLD_FRAC,
            baseline_name=r.GROIN_BASELINE_NAME)


# Section 12.4: everything the run writes into its folder besides the model state and figures
def write_run_outputs(r):
    # 12.4 Write

    _shoreface_depth_m = float(r.cascade._brie_coupler._brie.d_sf)

    # The shoreline matrix: what the figures, and a groin run's comparison, are drawn from
    np.save(write_path(r.RUN_DIR, "matrix", r.RUN_NAME), np.asarray(r.shoreline_m, dtype=float))

    # What the figures need from the model itself: overwash, and the years each road moved
    np.save(Path(r.RUN_DIR) / OVERWASH_NPY.format(run=r.RUN_NAME),
            np.array([np.asarray(_b3d.QowTS) for _b3d in r.cascade.barrier3d]).T)
    if getattr(r.cascade, "_roadways", None):
        _relocations = np.zeros((r.states, HATTERAS_DOMAINS.total_domains), dtype=bool)
        for _pad, _roadway in enumerate(r.cascade._roadways):
            # Not `getattr(...) or []`: the attribute is a numpy array
            _raw = getattr(_roadway, "_road_relocated_TS", None)
            _series = (np.asarray([], dtype=float) if _raw is None
                       else np.asarray(_raw, dtype=float))
            if _series.size:
                _n = min(_series.size, r.states)
                _relocations[:_n, _pad] = _series[:_n] > 0
        np.save(Path(r.RUN_DIR) / RELOCATIONS_NPY.format(run=r.RUN_NAME), _relocations)

    _rate_csv = str(write_path(r.RUN_DIR, "rate_csv", r.RUN_NAME))
    # Both estimators ship, with lrr_r2 beside them
    pd.DataFrame({
        "gis_domain": np.arange(HATTERAS_DOMAINS.first_gis_id,
                                HATTERAS_DOMAINS.last_gis_id + 1),
        "change_rate_m_yr": r.change_rate[r.real],
        "lrr_m_yr": r.model_lrr[r.real],
        "lrr_r2": r.model_lrr_r2[r.real],
    }).to_csv(_rate_csv, index=False)
    print(f"\nwrote                 {os.path.basename(_rate_csv)}")

    if r.ROAD_SUMMARY:
        _road_csv = str(write_path(r.RUN_DIR, "road_csv", r.RUN_NAME))
        pd.DataFrame(r.ROAD_SUMMARY).to_csv(_road_csv, index=False)
        print(f"                      {os.path.basename(_road_csv)}")

    # Skill against the section 8 target

    # Island-wide and interior spans both reported; SKILL is the LRR one, SKILL_ENDPOINT kept beside it
    SKILL = skill_vs_target(r.model_lrr, r.COASTSAT_TARGET, HATTERAS_DOMAINS)
    SKILL_ENDPOINT = skill_vs_target(r.change_rate, r.COASTSAT_TARGET,
                                     HATTERAS_DOMAINS)
    # Interior = GIS 2-89 against the surveyed target in every geometry
    for _skill, _rate in ((SKILL, r.model_lrr), (SKILL_ENDPOINT, r.change_rate)):
        _skill["mean_bias_reach_interior_m_yr"] = _skill["mean_bias_interior_m_yr"]
        _skill["rmse_reach_interior_m_yr"] = _skill["rmse_interior_m_yr"]
        _base = skill_vs_target(_rate, r.COASTSAT_TARGET_BASE, HATTERAS_DOMAINS,
                                interior_gis=SCORE_INTERIOR_GIS)
        for _key in ("mean_bias_interior_m_yr", "rmse_interior_m_yr",
                     "n_domains_interior"):
            _skill[_key] = _base[_key]
    print(f"\nSKILL vs CoastSat     model LRR - target LRR, m/yr")
    print(f"  island-wide         bias {SKILL['mean_bias_m_yr']:+.3f}   "
          f"RMSE {SKILL['rmse_m_yr']:.3f}   (n={SKILL['n_domains']})")
    print(f"  interior (GIS 2-89) bias {SKILL['mean_bias_interior_m_yr']:+.3f}   "
          f"RMSE {SKILL['rmse_interior_m_yr']:.3f}   "
          f"(n={SKILL['n_domains_interior']})")
    print(f"  endpoint estimator  bias "
          f"{SKILL_ENDPOINT['mean_bias_interior_m_yr']:+.3f}   "
          f"RMSE {SKILL_ENDPOINT['rmse_interior_m_yr']:.3f}   "
          f"(interior; the pre-LRR metric)")

    # Where a straight line summarises the trajectory poorly (nourished or barely moving domains)
    _lrr_r2_real = r.model_lrr_r2[r.real]
    _poor_fit = [int(_gis) for _gis, _r2 in zip(
        range(HATTERAS_DOMAINS.first_gis_id, HATTERAS_DOMAINS.last_gis_id + 1),
        _lrr_r2_real) if np.isfinite(_r2) and _r2 < r.LRR_R2_FLOOR]
    print(f"\nLRR FIT QUALITY       median r2 "
          f"{np.nanmedian(_lrr_r2_real):.3f}   "
          f"{len(_poor_fit)} of {_lrr_r2_real.size} domains below "
          f"{r.LRR_R2_FLOOR:.2f}")
    if _poor_fit:
        print(f"  a step, or no trend GIS "
              f"{', '.join(str(_g) for _g in _poor_fit[:12])}"
              f"{' ...' if len(_poor_fit) > 12 else ''}")

    # Run metadata: the scenario, plus what distinguishes this run

    # One structure renders both the .txt and the .json
    _GIT = git_provenance(REPO_DIR)
    _TIMESTAMP = timestamp()
    # The file that ran: the notebook sets r.GENERATED_BY; the script leaves the default
    GENERATED_BY = getattr(r, "GENERATED_BY", "HAT_hindcast.py")

    # What the source/sink preset actually was

    # The preset's name and its numbers are recorded separately
    _BE_NONZERO = sum(1 for _rate in r.DOMAIN_BE_RATES.values() if _rate)
    _BE_DIGEST = values_digest(r.DOMAIN_BE_RATES)
    _BE_EDGE_RATES = {f"rate_gis{_gis}_m_yr": r.DOMAIN_BE_RATES.get(_gis, 0.0)
                      for _gis in HATTERAS_BE_EDGE_DOMAINS}

    _META = {
        "identity": {
            "run_name": r.RUN_NAME,
            "timestamp": _TIMESTAMP,
            "generated_by": GENERATED_BY,
            "runtime_min": f"{r.RUN_SECONDS / 60:.1f}",
            # The three that drift across a batch: Cascade class, extractor version, commit
            "use_sandbox_cascade": (USE_SANDBOX_CASCADE,
                                    "cascade.cascade_groin, not cascade.cascade"),
            "cascade_patches": (list(r.CASCADE_PATCHES),
                                "obx-cascade functions replaced at run time (Hatteras/cascade_patch.py)"),
            # Product and version both: version numbers restart per product
            "topo_product": r.TOPO_PRODUCT,
            "topo_dune_version": r.TOPO_DUNE_VERSION,
            # The island offset version, which the run name does not carry
            "island_offset_version": r.ISLAND_OFFSET_VERSION,
            # The reach: "base" is GIS 1-90; an extended one is not a matrix run
            "geometry": (HATTERAS_GEOMETRY,
                         f"GIS {HATTERAS_DOMAINS.first_gis_id} to "
                         f"{HATTERAS_DOMAINS.last_gis_id}"),
            "run_kind": r.RUN_KIND,
            "run_tag": r.RUN_TAG,
            "parameter_file": (r.RUN_PARAMETER_FILE.name,
                               f"this run's copy of {r.PARAMETER_FILE}; CASCADE "
                               f"rewrites the copy, not the template"),
            "save_model_state": r.SAVE_MODEL_STATE,
            "git_commit": _GIT["commit"],
            "git_branch": _GIT["branch"],
            "git_dirty": (_GIT["dirty"],
                          "True: the commit alone does not reproduce this run"),
            # Barrier3D's branch is part of the model
            "barrier3d_branch": r.B3D["branch"],
            "barrier3d_commit": r.B3D["commit"],
            "barrier3d_dirty": r.B3D["dirty"],
            "barrier3d_route_overwash_fix": (r.B3D["route_overwash_fix"],
                                             "True: the loaded Barrier3D has the "
                                             "2026-09-24 route_overwash index fix"),
            "barrier3d_gap_momentum_fix": (r.B3D.get("gap_momentum_fix"),
                                           "True: DuneGaps, gap discharge slice and "
                                           "inundation momentum fixed (2026-09-28)"),
            "barrier3d_per_cell_ceiling": r.B3D.get("per_cell_ceiling"),
            # Read off the constructed model, not the template
            "dune_ceiling": ("per-cell, from the starting dunes"
                             if getattr(r.cascade.barrier3d[0], "_DuneCeilingFromStart", False)
                             else f"uniform Dmaxel {r.cascade.barrier3d[0].Dmaxel * 10 + r.MHW_ELEVATION:.2f} m NAVD88"),
            "storm_file": str(r.STORM_FILE.relative_to(r.HATTERAS_DATA_BASE)),
            # beach_dune_manager's 4 m bulldozer cap: what it clips (2026-09-28).
            "bdm_dune_cap_applies_to": getattr(_bdm_module, "DUNE_CAP_APPLIES_TO", "whole dune cell"),
        },
        "scenario": {label: value for label, value, _token in r.SCENARIO_SWITCHES},
        "period": {
            "start_year": r.START_YEAR,
            "end_year": r.END_YEAR,
            "first_model_year": r.FIRST_MODEL_YEAR,
            "last_model_year": r.LAST_MODEL_YEAR,
            "run_years": (r.RUN_YEARS, "transitions"),
            "annual_states": r.states,
        },
        "wave climate": {
            "wave_height_m": r.Hs,
            "wave_period_s": r.FIXED_WAVE_PERIOD,
            "wave_asymmetry": r.FIXED_WAVE_ASYMMETRY,
            "wave_angle_high_frac": r.FIXED_WAVE_ANGLE_HIGH_FRACTION,
            "shoreface_depth_m": (f"{_shoreface_depth_m:.2f}", "8.9 * Hs"),
        },
        "sea level": {
            "rslr_m_yr": r.SEA_LEVEL_RISE_RATE,
            "rslr_constant": r.SEA_LEVEL_CONSTANT,
        },
        "dunes and roadway": {
            "rmin / rmax": f"{r.RMIN[0]} / {r.RMAX[0]}",
            "growth_override": r.GROWTH_OVERRIDE_RAW or "none",
            "dune_design_ele_m_mhw": (r.DUNE_DESIGN_ELEVATION_M,
                                      "floored by roadway_manager"),
            "dune_min_ele_m_mhw": (r.DUNE_MINIMUM_ELEVATION_M,
                                   "floored by roadway_manager"),
            "road_width_m": r.ROAD_WIDTH,
            "road_dune_rebuild_trigger_m": (r.ROAD_DUNE_REBUILD_TRIGGER_M,
                                            "no road dune rebuild beyond this setback"),
            "relocations_enabled": r.ENABLE_HISTORICAL_ROAD_RELOCATIONS,
        },
        "source/sink": {
            "preset": r.SOURCE_SINK_PRESET,
            "background_erosion_on": (r.USE_BACKGROUND_EROSION,
                                      "implied by the preset, checked in 4.3"),
            "domains_specified": len(r.DOMAIN_BE_RATES),
            "nonzero_domains": _BE_NONZERO,
            "values_digest": (_BE_DIGEST,
                              "changes if any rate changed, interior included"),
            **_BE_EDGE_RATES,
        },
        "groin": {"enabled": r.GROIN_ENABLED},
        "skill": {
            "target": (f"CoastSat LOWESS {r.TARGET_WINDOW}-domain",
                       "raw means over GIS 1-"
                       f"{r.LOWESS_CONFIG.skip_southern_domains}"),
            "estimator": (r.RATE_ESTIMATOR,
                          "OLS slope through every annual state, matching the "
                          "target's own definition"),
            "mean_bias_m_yr": f"{SKILL['mean_bias_m_yr']:+.4f}",
            "rmse_m_yr": f"{SKILL['rmse_m_yr']:.4f}",
            "mean_bias_interior_m_yr": (
                f"{SKILL['mean_bias_interior_m_yr']:+.4f}",
                "GIS 2-89: the locked end domains excluded"),
            "rmse_interior_m_yr": f"{SKILL['rmse_interior_m_yr']:.4f}",
            "lrr_r2_median": (f"{np.nanmedian(_lrr_r2_real):.4f}",
                              "how well a line describes the modelled trajectory"),
            "lrr_r2_below_floor": (f"{len(_poor_fit)}",
                                   f"domains under r2 {r.LRR_R2_FLOOR:.2f}"),
            "endpoint_mean_bias_interior_m_yr": (
                f"{SKILL_ENDPOINT['mean_bias_interior_m_yr']:+.4f}",
                "the pre-LRR estimator; calibBE and groin M were fit on this"),
            "endpoint_rmse_interior_m_yr":
                f"{SKILL_ENDPOINT['rmse_interior_m_yr']:.4f}",
        },
        "verification": {
            "nourishment_ok": r.BN_REPORT["ok"],
            "roads_drowned": len(r.drowned),
            "roads_reloc_blocked": len(r.blocked),
        },
    }

    if r.cascade.groin is not None:
        _META["groin"].update({
            "kind": r.cascade.groin.kind,
            "strength_history": r.cascade.groin.strength_history,
            "spread_weights": list(r.cascade.groin.spread_weights),
            "mean_trapping_rate_m_yr": f"{r.cascade.groin.mean_trapping_rate_m_yr:.2f}",
            "updrift / downdrift": f"GIS {r.GROIN_UPDRIFT_GIS} / {r.GROIN_DOWNDRIFT_GIS}",
            "install_year": r.cascade.groin.install_year,
            "r_ipl_t0": f"{r.R_IPL:.4f}",
            "r_ipl_theta_deg": f"{r.R_IPL_THETA_DEG:.1f}",
            # The paired no-groin run, which the difference animation is drawn against
            "baseline_run": r.GROIN_BASELINE_NAME,
        })
        if r.GROIN_EXTENT is not None:
            _META["groin"]["measured_extent_m"] = (
                f"{r.GROIN_EXTENT['updrift_m']:.0f} updrift / "
                f"{r.GROIN_EXTENT['downdrift_m']:.0f} downdrift")

    # Cross-run row: the run's summary, which run_registry.scan_runs reads back for comparing the matrix
    _index_row = {
        "run_name": r.RUN_NAME,
        "timestamp": _TIMESTAMP,
        "start_year": r.START_YEAR,
        "end_year": r.END_YEAR,
        "source_sink_preset": r.SOURCE_SINK_PRESET,
        "be_nonzero_domains": _BE_NONZERO,
        "be_values_digest": _BE_DIGEST,
        **{f"be_{_key}": _value for _key, _value in _BE_EDGE_RATES.items()},
        "scenario": r.SCENARIO,
        "scenario_overridden": bool(r.SCENARIO_DEPARTURES),
        "groin_enabled": r.GROIN_ENABLED,
        "groin_strength": r.GROIN_STRENGTH if r.cascade.groin is not None else np.nan,
        "groin_end_year": r.GROIN_END_YEAR if r.cascade.groin is not None else np.nan,
        # The shoreline change the groin moved across its face, mean over its active years
        "groin_trapping_m_yr": (np.nan if r.cascade.groin is None
                                else r.cascade.groin.mean_trapping_rate_m_yr),
        "roadway_management": r.ENABLE_ROADWAY_MANAGEMENT,
        "relocations_enabled": r.ENABLE_HISTORICAL_ROAD_RELOCATIONS,
        "beach_dune_management": r.ENABLE_BEACH_DUNE_MANAGEMENT,
        "nourishment_fills": r.ENABLE_NOURISHMENT_FILLS,
        "bdm_domains": int(sum(r.BEACH_DUNE_MANAGEMENT_ON)),
        "nourishment_projects": len(r.BN_SCHEDULE_APPLIED.projects),
        "Hs_m": r.Hs,
        # Where the run is filed, spelled out
        "kind": r.RUN_KIND,
        "tag": r.RUN_TAG,
        "sandbags_on": r.ENABLE_SANDBAG_PLACEMENT,
        "rslr_m_yr": r.SEA_LEVEL_RISE_RATE,
        "annual_states": r.states,
        "run_complete": r.states == r.RUN_YEARS + 1,
        "rate_estimator": r.RATE_ESTIMATOR,
        "mean_bias_m_yr": SKILL["mean_bias_m_yr"],
        "rmse_m_yr": SKILL["rmse_m_yr"],
        "mean_bias_interior_m_yr": SKILL["mean_bias_interior_m_yr"],
        "rmse_interior_m_yr": SKILL["rmse_interior_m_yr"],
        "lrr_r2_median": float(np.nanmedian(_lrr_r2_real)),
        "lrr_r2_below_floor": len(_poor_fit),
        "endpoint_mean_bias_interior_m_yr":
            SKILL_ENDPOINT["mean_bias_interior_m_yr"],
        "endpoint_rmse_interior_m_yr": SKILL_ENDPOINT["rmse_interior_m_yr"],
        "nourishment_ok": r.BN_REPORT["ok"],
        "roads_drowned": len(r.drowned),
        "roads_reloc_blocked": len(r.blocked),
        "groin_extent_updrift_m": (r.GROIN_EXTENT["updrift_m"]
                                   if r.GROIN_EXTENT is not None else np.nan),
        "groin_extent_downdrift_m": (r.GROIN_EXTENT["downdrift_m"]
                                     if r.GROIN_EXTENT is not None else np.nan),
        "runtime_min": round(r.RUN_SECONDS / 60, 1),
        "use_sandbox_cascade": USE_SANDBOX_CASCADE,
        "topo_product": r.TOPO_PRODUCT,          # see the note at the json write
        "topo_dune_version": r.TOPO_DUNE_VERSION,
        "island_offset_version": r.ISLAND_OFFSET_VERSION,
        "geometry": HATTERAS_GEOMETRY,
        "rmse_reach_interior_m_yr": SKILL["rmse_reach_interior_m_yr"],
        "git_commit": _GIT["commit"][:12],
        "git_dirty": _GIT["dirty"],
        "barrier3d_commit": str(r.B3D["commit"])[:12],
        "barrier3d_gap_momentum_fix": r.B3D.get("gap_momentum_fix"),
        "dune_ceiling": ("per-cell" if getattr(r.cascade.barrier3d[0], "_DuneCeilingFromStart", False)
                         else "uniform"),
        "storm_file": r.STORM_FILE.name,
        "bdm_dune_cap": getattr(_bdm_module, "DUNE_CAP_APPLIES_TO", "whole dune cell"),
        "barrier3d_route_overwash_fix": r.B3D["route_overwash_fix"],
    }
    # The row goes into the metadata; scan_runs builds the cross-run table from every run's metadata
    _META[INDEX_SECTION] = _index_row
    _meta_txt, _meta_json = write_run_metadata(
        r.RUN_DIR, r.RUN_NAME, _META,
        header=[f"CASCADE run metadata -- generated by {GENERATED_BY}, "
                f"section 12.",
                "Companion .json holds the same values, machine-readable."])
    print(f"                      {_meta_txt.name}")
    print(f"                      {_meta_json.name}")
