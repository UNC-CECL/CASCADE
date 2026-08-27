#!/usr/bin/env python3
"""The printed reports the Hatteras hindcast emits at each section.

WHY THIS MODULE EXISTS
    The reports are the run's evidence: they are what says a switch reached
    the module it names, what the sediment budget implies, and where the
    model disagrees with the survey. They are not debug output, so none of
    them were dropped -- but ~400 lines of f-strings sat duplicated between
    `HAT_hindcast_1984_2024.ipynb` and its headless mirror
    `HAT_hindcast_1984_2024.py`, and every wording change had to be made
    twice or the two runs stopped saying the same thing.

    Every function here prints exactly what the inline block printed. The
    values still come from the calling file -- these take arguments, never
    module globals -- so the notebook remains the place that decides what is
    reported, and this is only where the formatting lives.

ONE RULE
    A report never computes a result the run depends on. `run_units_check`
    is the single exception and it is named for it: it loops the domains and
    returns whether anything failed, because the loop existed only to be
    printed. Everything else takes already-built objects.

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import numpy as np

from cascade_pipeline.hindcast import (
    check_domain_units,
    groin_differential,
    groin_trapping_schedule,
    implied_interception_m3_yr,
)

__all__ = [
    "path_inventory", "run_units_check", "scenario_report", "storm_report",
    "background_erosion_report", "roadway_report", "road_audit_report",
    "beach_dune_report", "groin_report", "scenario_summary_report",
    "coastsat_report", "figure_config_report", "pre_run_report",
    "run_length_report", "nourishment_report", "frozen_setbacks_report",
    "roadway_outcome_report", "groin_extent_report", "target_misfit_report",
]


def path_inventory(pairs):
    """Prints name/exists/path for each (label, Path) pair."""
    for name, path in pairs:
        print(f"{name:<20} {'ok' if path.exists() else 'MISSING':<8} {path}")


# =============================================================================
# SECTION 2 -- UNITS
# =============================================================================

def run_units_check(elevation_paths, dune_paths, contract, geometry):
    """Runs the units contract over every real domain and reports the tally.

    Args:
        elevation_paths, dune_paths: Padded lists of .npy paths.
        contract: Mapping from load_barrier3d_contract.
        geometry: DomainGeometry.

    Returns:
        True if every check passed on every domain.
    """
    results = {}
    for gis_id in range(geometry.first_gis_id, geometry.last_gis_id + 1):
        pad = geometry.gis_to_pad(gis_id)
        checks = check_domain_units(np.load(elevation_paths[pad]),
                                    np.load(dune_paths[pad]), contract)
        for name, passed in checks.items():
            results.setdefault(name, []).append((gis_id, passed))

    print(f"{'Check':<46} {'domains passing':>15}")
    failed_any = False
    for name, outcomes in results.items():
        failures = [gis_id for gis_id, passed in outcomes if not passed]
        failed_any = failed_any or bool(failures)
        print(f"  {name:<44} {len(outcomes) - len(failures):>6}"
              f"/{len(outcomes)}"
              + (f"   FAILED: {failures[:8]}" if failures else ""))

    print("\n" + ("SOME CHECKS FAILED - do not run the model on these arrays"
                  if failed_any else
                  "All units and shapes match what Barrier3D expects."))
    return not failed_any


# =============================================================================
# SECTION 3 -- SCENARIO
# =============================================================================

def scenario_report(*, scenario, scenarios, departures, roadway_on,
                    relocations_on, relocations_forced_off, beach_dune_on,
                    fills_on, fills_forced_off, period_expects_nourishment,
                    groin_enabled, run_name_preview, input_files):
    """Prints the scenario, what it expands to, and the predicted run name.

    Args:
        scenario: Selected scenario key.
        scenarios: The full scenario preset table.
        departures: Mapping of key -> (wanted, got) for overridden switches.
        roadway_on, relocations_on, beach_dune_on, fills_on, groin_enabled:
            The resolved switches.
        relocations_forced_off, fills_forced_off: Whether the switch was
            forced off by its parent module rather than chosen.
        period_expects_nourishment: PERIOD["enable_nourishment"], printed
            beside the switch that does reach the model.
        run_name_preview: The name predicted from the switches.
        input_files: Sequence of (label, Path) for the period's inputs.
    """
    print(f"\nSCENARIO              {scenario!r}"
          + (f"   OVERRIDDEN: {', '.join(sorted(departures))}"
             if departures else ""))
    for key, preset in scenarios.items():
        print(f"  {'->' if key == scenario else '  '} {key:<18}"
              f"road={str(preset['roadway']):<5} "
              f"bdm={str(preset['beach_dune']):<5} "
              f"fills={str(preset['fills'])}")
    for key, (want, got) in sorted(departures.items()):
        print(f"  override            {key}: scenario says {want}, "
              f"this run uses {got}")

    print("\nMANAGEMENT            what the scenario above expands to")
    print(f"  roadway_manager     {'on' if roadway_on else 'OFF'}")
    print(f"    relocations       {'on' if relocations_on else 'off'}"
          + ("   FORCED off: roadway_manager is off"
             if relocations_forced_off else ""))
    print(f"  beach_dune_manager  {'on' if beach_dune_on else 'OFF'}")
    print(f"    nourishment fills {'on' if fills_on else 'off'}"
          + ("   FORCED off: beach_dune_manager is off"
             if fills_forced_off else ""))
    # PERIOD["enable_nourishment"] describes the period; nothing reads it. It
    # is printed beside the switch that does reach the model so the two cannot
    # quietly disagree -- the volume there is a legacy default that section 6
    # overwrites. Flagged in one direction only: a period that expects fill and
    # is not getting it is a scenario choice worth restating, whereas fills
    # left on in a period with no projects withholds nothing.
    print(f"    period expects    {period_expects_nourishment}"
          + ("   <- the period has fill and this run suppresses it"
             if period_expects_nourishment and not fills_on else ""))
    print(f"  groin               {'on' if groin_enabled else 'off'}"
          f"   (paired baseline: run each scenario with False first)")

    print(f"\nRUN_NAME (predicted)  {run_name_preview!r}")
    print(f"                      confirmed against the built modules in 7.5")
    print()
    for name, path in input_files:
        print(f"  {name:<14} {'ok' if path.exists() else 'MISSING':<8} "
              f"{path.name}")


# =============================================================================
# SECTION 4 -- FORCINGS
# =============================================================================

def storm_report(*, storms, storm_file, run_years):
    """Prints the storm series extent and how much of it the run reaches."""
    print(f"\n{storm_file.name}: {len(storms)} storms")
    print(f"  time steps    {storms['time'].min():.0f} to "
          f"{storms['time'].max():.0f}  (run is {run_years} years)")
    print(f"  Rhigh         {storms['Rhigh_m'].min():.2f} to "
          f"{storms['Rhigh_m'].max():.2f} m")
    print(f"  duration      {storms['duration'].min():.0f} to "
          f"{storms['duration'].max():.0f} hours")
    print(f"  storms/year   {len(storms) / run_years:.1f} mean")

    beyond_run = (storms["time"] > run_years + 1).sum()
    if beyond_run:
        print(f"  NOTE: {beyond_run} storms land past model year "
              f"{run_years + 1}")


def background_erosion_report(*, preset, start_year, domain_rates, rates,
                              use_background_erosion, geometry,
                              rates_2004_are_placeholder):
    """Prints the source/sink preset, its coverage, and its range."""
    print(f"\npreset {preset!r}, period {start_year}: "
          f"{len(domain_rates)} domains specified, "
          f"{sum(1 for r in rates if r != 0)} non-zero of "
          f"{geometry.total_domains}")
    if domain_rates:
        print(f"  rate range     {min(domain_rates.values()):+.1f} to "
              f"{max(domain_rates.values()):+.1f} m/yr")
    else:
        print("  rate range     no domains specified -- 0.0 m/yr everywhere")
    print(f"  USE_BACKGROUND_EROSION = {use_background_erosion}")
    if (start_year == 2004 and preset == "calibBE"
            and rates_2004_are_placeholder):
        print("  WARNING: 2004 calibrated rates are a copy of the 1984 fit")


# =============================================================================
# SECTION 5 -- ROADWAY
# =============================================================================

def roadway_report(*, setback_file, setbacks, missing_setbacks,
                   elevation_file, elevations, missing_elevations,
                   road_slice, config, roadway_on, management_on,
                   first_road_gis, last_road_gis, community_zones,
                   road_events, relocations_enabled):
    """Prints the roadway forcing, the managed footprint, and the events."""
    print(f"\nsetback file    {setback_file.name}")
    print(f"  {setbacks[road_slice].min():.0f}"
          f"-{setbacks[road_slice].max():.0f} m"
          + (f"   MISSING {missing_setbacks}" if missing_setbacks else ""))
    print(f"elevation file  {elevation_file.name}  (period-independent)")
    print(f"  {elevations[road_slice].min():.2f}"
          f"-{elevations[road_slice].max():.2f} m MHW"
          f"   buffers {config.elevation_fallback_m} m"
          + (f"   MISSING {missing_elevations}" if missing_elevations else ""))
    if roadway_on:
        print(f"managed         {sum(management_on)} domains "
              f"({last_road_gis - first_road_gis + 1} "
              f"carry road, minus community zones "
              f"{list(community_zones)})")
    else:
        print(f"managed         0 domains -- ENABLE_ROADWAY_MANAGEMENT is "
              f"False (section 3)")
        print(f"                no bulldozing, no dune rebuild, no setback "
              f"tracking, no road drowning")
    for event in road_events:
        print(f"event {event.year}      {type(event).__name__:<15} "
              f"{event.note}")
    print(f"ENABLE_HISTORICAL_ROAD_RELOCATIONS = {relocations_enabled}")


def road_audit_report(*, audit, summary):
    """Prints which roads would not survive year one, and why.

    A domain whose road starts in water is not a warning: roadway_manager
    sets _drown_break and returns immediately on every later year, so the
    domain becomes an unmanaged barrier wearing a road label.
    """
    print(f"\n{summary['n_domains']} road domains, "
          f"{summary['n_managed']} managed\n")
    if summary["blocking"]:
        print("BLOCKING - these would corrupt the run rather than drown it:")
        for gis, wall in summary["blocking"]:
            print(f"  GIS {gis}: {wall}")
        print()
    print(f"drowns at t=0, managed   : {summary['drowning']}")
    print(f"drowns at t=0, unmanaged : {summary['drowning_unmanaged']}"
          f"   (community zones; no manager runs there)\n")

    bad = [r for r in audit if r["drowns"] and r["managed"]]
    if bad:
        print(f"{'GIS':>4} {'setback':>9} {'road rows':>11} {'%sea':>6} "
              f"{'%bay':>6} {'%road wet':>10} {'land ends min/med':>18}")
        for r in bad:
            rows = f"{r['road_start']}-{r['road_end'] - 1}"
            ends = f"{r['width_min']} / {r['width_median']}"
            print(f"{r['gis']:>4} {r['setback_m']:>8.0f}m {rows:>11} "
                  f"{r['sea_water'] * 100:>5.0f}% {r['bay_water'] * 100:>5.0f}% "
                  f"{r['road_cells_water'] * 100:>9.0f}% {ends:>18}")


# =============================================================================
# SECTION 6 -- BEACH/DUNE MANAGER
# =============================================================================

def beach_dune_report(*, start_year, end_year, beach_dune_enabled,
                      fills_enabled, roadway_enabled, config, management_on,
                      overwash_to_dune, community_zones, schedule,
                      schedule_applied, audit, double_managed, geometry):
    """Prints the module footprint, the fill schedule, and the audit."""
    print(f"\nperiod                {start_year}-{end_year}")
    if beach_dune_enabled:
        print(f"overwash_filter       "
              f"{config.community_overwash_filter_pct:.0f}% in "
              f"community zones, "
              f"{config.default_overwash_filter_pct:.0f}% "
              f"elsewhere   (PERCENT, not a fraction)")
        print(f"overwash_to_dune      {overwash_to_dune:.0f}%")
        print(f"module on             {sum(management_on)} domains "
              f"= community zones {list(community_zones)} "
              f"UNION nourished {list(schedule.nourished_gis)}")
        if not fills_enabled and schedule.projects:
            print(f"                      the footprint still carries the "
                  f"project extents: the fills are off, the module is not")
    else:
        print(f"module on             0 domains -- "
              f"ENABLE_BEACH_DUNE_MANAGEMENT is False (section 3)")
        print(f"                      no overwash filtering, no fixed dune "
              f"line, no width-drowning check, no fill")
        print(f"                      overwash_filter handed to CASCADE is "
              f"all zero, so the array says what the run does")

    print(f"\nnourishment schedule  "
          f"{len(schedule_applied.projects)} project(s) applied"
          + (f"   {len(schedule.projects)} in period SUPPRESSED by "
             f"ENABLE_NOURISHMENT_FILLS=False"
             if not fills_enabled and schedule.projects else ""))
    for project in schedule_applied.projects:
        per_m = project.volume_m3_per_m(geometry.domain_spacing_m)
        print(f"  {project.year}  {project.name:<26} "
              f"GIS {project.gis_domains[0]}-{project.gis_domains[-1]:<3} "
              f"{project.volume_cubic_yards:>9,.0f} cy -> "
              f"{per_m:>7.1f} m^3/m per domain   "
              f"(TS index {schedule_applied.time_index(project.year)})")
    for project, reason in schedule.skipped:
        print(f"  --    {project.name:<26} skipped: {reason}")

    print(f"\npre-flight audit      {'ok' if audit['ok'] else 'PROBLEMS'}")
    if audit["module_off"]:
        print(f"  scheduled but module OFF (fill would be dropped): "
              f"{audit['module_off']}")
    for row in audit["implausible_volume"]:
        print(f"  implausible volume: GIS {row['gis']} {row['year']} "
              f"{row['volume_m3_per_m']:.1f} m^3/m")
    for note in audit["notes"]:
        print(f"  note: {note}")

    if double_managed:
        print(f"\nBOTH MODULES ON       GIS {list(double_managed)} "
              f"({len(double_managed)} domains)")
        print("  Overwash is removed twice over: RoadwayManager bulldozes and")
        print("  rewrites DomainTS, then BeachDuneManager filters what "
              "survived.")
        print("  BeachDuneManager also holds dune_migration_on False, so")
        print("  ShorelineChangeTS stays 0 -- the value RoadwayManager reads "
              "to")
        print("  decrement the setback. NC-12 stops retreating here.")
        print("  Run as published; verified against the finished run in "
              "section 12.")
    elif not (roadway_enabled and beach_dune_enabled):
        print("\nBOTH MODULES ON       none - at least one module is off, so "
              "the overlap cannot arise")
    else:
        print("\nBOTH MODULES ON       none - the two footprints are disjoint "
              "this period")


# =============================================================================
# SECTION 7 -- GROIN
# =============================================================================

def groin_report(*, enabled, callback, updrift_gis, downdrift_gis,
                 install_year, start_year, end_year, geometry,
                 trapping_rate_m_yr, m_provenance, deterioration_mode,
                 deterioration_delay_years, deterioration_ramp_years,
                 deterioration_fraction, profile_height_candidates_m,
                 reach_transport_loss_m3_yr, reach_transport_citation,
                 reach_transport_caveat, source_sink_preset, domain_be_rates):
    """Prints the groin configuration, its sediment budget, and the overlap.

    The double-count note is reported rather than corrected: the calibrated
    background-erosion rates were fit to CoastSat spanning the functional-groin
    era, so the groin's signature is plausibly already in them.
    """
    up_pad = geometry.gis_to_pad(updrift_gis)
    down_pad = geometry.gis_to_pad(downdrift_gis)
    years, m_eff = groin_trapping_schedule(callback, start_year, end_year)

    print(f"\ngroin                 {'ENABLED' if enabled else 'disabled'}"
          f"   (module cascade.groin, hook cascade_groin.py pre-AST)")
    print(f"structure             updrift GIS {updrift_gis} (pad {up_pad}) "
          f"/ downdrift GIS {downdrift_gis} (pad {down_pad}), "
          f"installed {install_year}")
    print(f"trapping rate M       {trapping_rate_m_yr:.0f} m/yr")
    print(f"  provenance          {m_provenance}")
    print(f"deterioration         {deterioration_mode}, "
          f"onset {callback.deterioration_year:.0f} "
          f"(+{deterioration_delay_years:.0f} yr), "
          f"ramp {deterioration_ramp_years:.0f} yr, "
          f"floor {deterioration_fraction:.2f}")
    print(f"  M over {start_year}-{end_year}    "
          f"{m_eff[0]:.1f} -> {m_eff[-1]:.1f} m/yr   "
          f"(mean {m_eff.mean():.1f})")
    if install_year <= start_year:
        print(f"  note                install {install_year} predates the "
              f"run: active every step, and {start_year - install_year} yr "
              f"of fillet is already in the initial shoreline")
    if deterioration_fraction >= 0.5:
        print(f"  TENSION             floor {deterioration_fraction:.2f} "
              f"leaves the structure at "
              f"{deterioration_fraction * 100:.0f}% strength after failure; "
              f"the 1996/2003 record does not support that")

    print(f"\nSEDIMENT BUDGET       reference {reach_transport_loss_m3_yr:,.0f} "
          f"m3/yr")
    print(f"  source              {reach_transport_citation}")
    print(f"  CAVEAT              {reach_transport_caveat}")
    print(f"  implied transfer at M = {trapping_rate_m_yr:.0f} m/yr over "
          f"dy = {geometry.domain_spacing_m:.0f} m:")
    for h in profile_height_candidates_m:
        vol = implied_interception_m3_yr(trapping_rate_m_yr, h, geometry)
        pct = 100.0 * vol / reach_transport_loss_m3_yr
        print(f"    profile height {h:5.1f} m ->  {vol:>9,.0f} m3/yr "
              f"= {pct:5.1f}% of the reach budget"
              + ("   BREACH" if pct > 50 else ""))
    print("  Not corrected. Section 11 resolves the profile height from the "
          "constructed model")
    print("  rather than the yaml, and section 12 reports the resulting "
          "misfit.")

    print(f"\nSOURCE/SINK OVERLAP   preset {source_sink_preset!r}, "
          f"period {start_year}   ((-) erosion / (+) accretion)")
    for label, gis in (("updrift  ", updrift_gis),
                       ("downdrift", downdrift_gis)):
        be = domain_be_rates.get(gis, 0.0)
        dipole = (-trapping_rate_m_yr if gis == updrift_gis
                  else +trapping_rate_m_yr)
        print(f"  {label} GIS {gis:<3} BE {be:+6.2f} m/yr    "
              f"groin x_s_dt {dipole:+7.1f} m/yr")
    print("  The calibrated BE rates were fit to observed CoastSat LRR "
          "spanning the")
    print("  functional-groin era, so the groin's signature is plausibly "
          "counted twice.")
    print("  Reported pending the source/sink re-analysis; not corrected "
          "here.")


def scenario_summary_report(*, scenario, departures, switches, run_name_base,
                            double_managed, groin_enabled, updrift_gis,
                            beach_dune_on_updrift):
    """Prints every switch, the derived run name, and the two known overlaps.

    Args:
        switches: The (label, value, token) list the name was derived from.
        beach_dune_on_updrift: Whether beach_dune_manager runs on the updrift
            domain. Keyed on the module flag, not the schedule: dune_migration
            is held False wherever the manager runs, fill year or not.
    """
    print(f"\nSCENARIO              {scenario!r}"
          + (f"   OVERRIDDEN: {', '.join(sorted(departures))}"
             if departures else ""))
    for label, value, token in switches:
        print(f"  {label:<24} {str(value):<22} "
              f"{'-> ' + token if token else ''}")
    print(f"\nRUN_NAME_BASE         {run_name_base!r}")
    print(f"                      matches the section 3 preview (asserted)")
    if double_managed:
        print(f"  both managers on    GIS {list(double_managed)} "
              f"-- see section 6")
    if groin_enabled and beach_dune_on_updrift:
        print(f"  groin + beach/dune  updrift GIS {updrift_gis} is inside "
              f"the beach_dune_manager footprint;")
        print(f"                      beach_dune_manager holds "
              f"dune_migration_on False there")


# =============================================================================
# SECTION 8 -- COASTSAT TARGET
# =============================================================================

def coastsat_report(*, target, active, target_window, loess_config, geometry,
                    cs_series, updrift_gis, downdrift_gis):
    """Prints the target curve, the unsmoothed southern zone, and the dipole."""
    n_loess = int((target["source"].str.startswith("LOESS")).sum())
    n_raw = len(target) - n_loess
    missing = target["target_lrr_m_yr"].isna().sum()

    print(f"\nactive dataset        {active['label']}")
    print(f"target window         LOESS {target_window}-domain "
          f"({target_window * geometry.domain_spacing_m / 1000:.1f} km)"
          f"   -- rate_comparison uses max(window_domains), asserted above")
    print(f"COASTSAT_TARGET       {len(target)} domains: "
          f"{n_raw} raw mean (D1-{loess_config.skip_southern_domains}), "
          f"{n_loess} LOESS"
          + (f", {missing} with no transects" if missing else ""))
    print(f"  range               "
          f"{target['target_lrr_m_yr'].min():+.2f} to "
          f"{target['target_lrr_m_yr'].max():+.2f} m/yr")
    print("  This is the curve the 'calibBE' source/sink preset was fit "
          "against (section 4.3).")

    print(f"\nUNSMOOTHED ZONE       D1-{loess_config.skip_southern_domains}, "
          f"raw per-domain means -- no LOESS line here")
    for row in target[target["gis_domain"]
                      <= loess_config.skip_southern_domains].itertuples():
        mark = ""
        if row.gis_domain == updrift_gis:
            mark = "  <- groin updrift"
        elif row.gis_domain == downdrift_gis:
            mark = "  <- groin downdrift"
        print(f"    D{row.gis_domain:<3} {row.target_lrr_m_yr:+7.2f} m/yr   "
              f"n={row.n_transects:<3}{mark}")

    print(f"\nGROIN DIFFERENTIAL    observed updrift D{updrift_gis} minus "
          f"downdrift D{downdrift_gis}")
    print("  positive = updrift doing better than downdrift = groin-like")
    for cs in cs_series:
        d = groin_differential(cs, updrift_gis, downdrift_gis)
        verdict = ("supports the dipole direction" if d["differential"] > 0
                   else "OPPOSITE to the imposed dipole")
        print(f"  {d['label']:<28} {'(active)' if d['active'] else '        '} "
              f"up {d['updrift']:+6.2f}  down {d['downdrift']:+6.2f}  "
              f"diff {d['differential']:+6.2f} m/yr")
        print(f"    {'':<28}          n = {d['n_updrift']} / "
              f"{d['n_downdrift']} transects   {verdict}")


# =============================================================================
# SECTION 9 -- FIGURE CONFIGURATION
# =============================================================================

def figure_config_report(*, annotations, loess_config, gif_config, gif_jobs,
                         flip_sign_model, real_domains_only, groin_enabled,
                         baseline_npy, baseline_name, output_base_dir):
    """Prints the figure configuration and whether the difference GIFs can run."""
    print(f"\nannotations           "
          f"{len(annotations.town_spans)} towns, "
          f"{len(annotations.piers)} piers, "
          f"{len(annotations.groins)} groin(s), "
          f"{len(annotations.shoal_zones)} shoal zones")
    print(f"loess_config          windows {loess_config.window_domains}, "
          f"skip D1-{loess_config.skip_southern_domains}   (section 8's)")
    print(f"gif_config            {gif_config.fps} fps, stride "
          f"{gif_config.year_stride}, ocean at "
          f"{'bottom' if gif_config.ocean_at_bottom else 'top'}, "
          f"save_matrix={gif_config.save_matrix}")
    print(f"sign convention       flip_sign_model={flip_sign_model} "
          f"(up = seaward), real domains only={real_domains_only}")

    print(f"\nGIF_JOBS              {len(gif_jobs)} job(s)")
    diff_jobs = [j for j in gif_jobs if j.get("mode") == "difference"]
    for job in gif_jobs:
        print(f"  {job['range']:<11} {job['mode']:<13}"
              + (f" pad={job['pad']}" if "pad" in job else ""))

    print("\nbaseline              ", end="")
    if not groin_enabled:
        print("n/a - this run has the groin off, so it IS the baseline")
        if diff_jobs:
            print(f"  {len(diff_jobs)} difference job(s) will be skipped by "
                  f"make_all_shoreline_gifs")
    elif baseline_npy:
        print(f"found\n  {baseline_name}\n  {baseline_npy}")
    else:
        print(f"NOT FOUND\n  looked for {baseline_name}")
        print(f"  under {output_base_dir}")
        if diff_jobs:
            print(f"  {len(diff_jobs)} difference job(s) will be SKIPPED. "
                  f"Run once with GROIN_ENABLED=False to create it.")


# =============================================================================
# SECTION 11 -- PRE-RUN
# =============================================================================

def pre_run_report(*, run_name, run_dir, run_years, start_year, end_year,
                   geometry, roadway_on, beach_dune_on, wave_height,
                   d_sf_m, h_b_m, profile_height_m, berm_floor_m,
                   dune_design_elevation_m, dune_minimum_elevation_m,
                   groin_callback, groin_enabled, r_ipl,
                   predicted_amplitude_m, predicted_extent_domains,
                   predicted_extent_m, reach_transport_loss_m3_yr):
    """Prints what was built, the resolved profile, and the fillet prediction.

    The prediction is printed BEFORE the run because a prediction printed
    after one is not a prediction. Amplitude was tuned; extent was not.
    """
    print(f"\nbuilt                 {run_name}")
    print(f"  run_dir             {run_dir}")
    print(f"  {run_years} transitions -> {run_years + 1} annual states "
          f"({start_year}-{end_year})")
    print(f"  domains             {geometry.total_domains} padded, "
          f"{sum(roadway_on)} roadway, "
          f"{sum(beach_dune_on)} beach/dune")

    print(f"\nSHOREFACE             Hs {wave_height} m -> d_sf = 8.9*Hs = "
          f"{d_sf_m:.2f} m (brie.py:270)")
    print(f"  h_b at t=0          {h_b_m:.2f} m")
    print(f"  active profile      {profile_height_m:.2f} m  "
          f"(resolves the section 7 volume diagnostic)")

    print(f"\nDUNE THRESHOLDS       passed -> floor roadway_manager will apply")
    print(f"  design              {dune_design_elevation_m:.2f} -> "
          f"{max(dune_design_elevation_m, berm_floor_m + 1.0):.2f} m MHW")
    print(f"  minimum             {dune_minimum_elevation_m:.2f} -> "
          f"{max(dune_minimum_elevation_m, berm_floor_m + 0.3):.2f} m MHW")
    print(f"  BermEl * 10         {berm_floor_m:.2f} m   "
          + ("(consistent with the metre reading)" if berm_floor_m < 5
             else "SUSPECT -- see the notebook markdown on BermEl units"))

    if groin_callback is not None:
        vol = implied_interception_m3_yr(groin_callback.M, profile_height_m,
                                         geometry)
        pct = 100.0 * vol / reach_transport_loss_m3_yr
        print(f"\nGROIN                 attached, M = {groin_callback.M:.0f} "
              f"m/yr, updrift pad {groin_callback.updrift_pad} / "
              f"downdrift pad {groin_callback.downdrift_pad}")
        print(f"  r_ipl (t=0)         {r_ipl:.4f}  (shore-normal, "
              f"brie.py:1294)")
        print(f"  PREDICTED amplitude {predicted_amplitude_m:.1f} m       "
              f"= M / (4 * r_ipl)")
        print(f"  PREDICTED extent    {predicted_extent_domains:.1f} domains "
              f"({predicted_extent_m:.0f} m)   = dy * sqrt(2 * r_ipl * t)")
        print(f"  Extent does not depend on M. Section 12 checks the emergent")
        print(f"  extent against this number, written down before the run.")
        print(f"  implied transfer    {vol:,.0f} m3/yr = {pct:.0f}% of the "
              f"{reach_transport_loss_m3_yr:,.0f} m3/yr reach budget")
        if pct > 50:
            print(f"                      BREACH, reported not corrected "
                  f"(section 7)")
    else:
        print(f"\nGROIN                 not attached (GROIN_ENABLED = "
              f"{groin_enabled})")
        print(f"  r_ipl (t=0)         {r_ipl:.4f}")


# =============================================================================
# SECTION 12 -- VERIFY
# =============================================================================

def run_length_report(*, states, run_years):
    """Prints how many annual states the run produced."""
    print(f"\nRUN LENGTH            {states} annual states for {run_years} "
          f"transitions")
    if states != run_years + 1:
        print(f"  INCOMPLETE          expected {run_years + 1}; the run "
              f"stopped early (b3d_break or drowning)")


def nourishment_report(*, report, run_years):
    """Prints the model's own record of what fill it received.

    Read from each manager's _nourishment_volume_TS, not from the schedule's
    intent -- the failure this catches produced plausible-looking output for
    a long time.
    """
    print(f"  denominator         {run_years} yr == states - 1   (asserted)")
    print(f"\nNOURISHMENT           {'OK' if report['ok'] else 'FAILED'}   "
          f"(read from each manager's _nourishment_volume_TS)")
    for row in report["confirmed"]:
        print(f"  ok       GIS {row['gis']:>3} {row['year']}  "
              f"{row['actual_m3_per_m']:>7.1f} m^3/m")
    for row in report["wrong_volume"]:
        print(f"  WRONG    GIS {row['gis']:>3} {row['year']}  expected "
              f"{row['expected_m3_per_m']:.1f}, got "
              f"{row['actual_m3_per_m']:.1f} m^3/m")
    for row in report["missing"]:
        print(f"  MISSING  GIS {row['gis']:>3} {row['year']}  expected "
              f"{row['expected_m3_per_m']:.1f} m^3/m, never applied")
    for row in report["unexpected"]:
        print(f"  EXTRA    GIS {row['gis']:>3} {row['year']}  "
              f"{row['volume_m3_per_m']:.1f} m^3/m, not scheduled")
    for note in report["notes"]:
        print(f"  note: {note}")


def frozen_setbacks_report(*, double_managed, rows):
    """Prints whether the setbacks section 6 predicted would freeze did."""
    if not double_managed:
        return
    print(f"\nFROZEN SETBACKS       predicted in section 6 for GIS "
          f"{list(double_managed)}")
    for row in rows:
        print(f"  GIS {row['gis']:>3}  {row['setback_start_m']:.0f} -> "
              f"{row['setback_end_m']:.0f} m"
              f"{'   FROZEN as predicted' if row['frozen'] else ''}")


def roadway_outcome_report(*, summary):
    """Prints which managed roads drowned and which relocations were blocked.

    Returns:
        The (drowned, relocation_blocked) row lists. Section 12.4 counts both
        into run_index.csv, so they are returned rather than recomputed.
    """
    drowned = [r for r in summary if r["drowned"]]
    blocked = [r for r in summary if r["relocation_blocked"]]
    print(f"\nROADWAY               {len(summary)} managed domains, "
          f"{len(drowned)} drowned, {len(blocked)} relocation-blocked")
    for row in drowned[:10]:
        print(f"  drowned  GIS {row['gis']:>3}  last managed "
              f"{row['last_managed_year']}  {row['reason']}")
    if len(drowned) > 10:
        print(f"  ... and {len(drowned) - 10} more (full list in the CSV)")
    return drowned, blocked


def groin_extent_report(*, extent, threshold_frac, baseline_name,
                        predicted_extent_domains, predicted_extent_m):
    """Prints the emergent fillet extent against the pre-registered prediction."""
    print(f"\nGROIN EXTENT          measured against {baseline_name}")
    print(f"  peak effect         {extent['peak_m']:.1f} m   "
          f"(threshold {threshold_frac:.0%} = "
          f"{extent['threshold_m']:.1f} m)")
    print(f"  updrift             {extent['updrift_domains']} domains "
          f"({extent['updrift_m']:.0f} m)")
    print(f"  downdrift           {extent['downdrift_domains']} domains "
          f"({extent['downdrift_m']:.0f} m)")
    print(f"  PREDICTED (sec 11)  "
          f"{predicted_extent_domains:.1f} domains "
          f"({predicted_extent_m:.0f} m)")
    print(f"  Amplitude was tuned; extent was not. This is the part that "
          f"could have failed.")


def target_misfit_report(*, target_m, observed_change_m, shoreline_m,
                         end_year, geometry, raw_offset_dir):
    """Prints the model's misfit against the surveyed end-year dune line."""
    print()
    if target_m is None:
        print(f"target                none -- no {end_year} dune-line survey "
              f"in {raw_offset_dir.name}/; GIFs draw the model alone")
        return
    real = slice(geometry.start_real_index, geometry.end_real_index)
    modeled = (shoreline_m[-1] - shoreline_m[0])[real]
    misfit = modeled - observed_change_m
    print(f"target                observed {end_year} dune line "
          f"({np.count_nonzero(~np.isnan(observed_change_m))} real domains)")
    print(f"  observed change     mean {np.nanmean(observed_change_m):+7.1f} m   "
          f"({np.nanmin(observed_change_m):+.1f} to "
          f"{np.nanmax(observed_change_m):+.1f})   + = landward")
    print(f"  modeled change      mean {modeled.mean():+7.1f} m   "
          f"({modeled.min():+.1f} to {modeled.max():+.1f})")
    print(f"  misfit              mean {np.nanmean(misfit):+7.1f} m   "
          f"RMSE {np.sqrt(np.nanmean(misfit ** 2)):.1f} m")
