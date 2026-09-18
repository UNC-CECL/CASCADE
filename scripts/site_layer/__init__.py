"""
site_layer -- what and where Hatteras is.

The seam between the site-agnostic `cascade_pipeline/` package and the task
folders under `scripts/`, which all consume Hatteras' own paths, place names
and house style. Six modules, each answering one question, each answering it
in one place:

    hatteras_site_config     what is this place -- domain geometry, town
                             spans, periods, BE presets, road events,
                             nourishment projects
    hat_topo_version         which Barrier3D domains -- product and version
    hat_elevation_products   which elevation product and stage
    hat_extension_domains    which alongshore reach, and the 500 m bins
                             beyond the 90 surveyed domains
    hat_observed_rates       where the observed shoreline is, and the rate
                             fits the model is graded against
    hat_figure_style         what a Hatteras figure looks like

WHY THE INDIRECTION EXISTS
    A hand-built path fails silently. Four road scripts once hardcoded
    `2009-dune-topo/2009_v3`; when the dune windows were re-picked into `v4`
    they kept reading `v3` interiors while consuming `v4` setbacks, and 18
    domains -- two of the three managed roadways among them -- had their drown
    verdicts computed on the wrong grid, with nothing raised. Each location is
    resolved once, here, and a name that is not on disk is an immediate, loud
    error listing what is.

    So never rebuild one of these paths by hand. Call `topo_dirs()`,
    `domain_arrays()`, `array_path()` or `product()` and let it raise.

THIS PACKAGE IS DELIBERATELY NOT RE-EXPORTED
    Import the module you need, not the package:

        from site_layer.hat_topo_version import topo_dirs
        from site_layer.hatteras_site_config import HATTERAS_DOMAINS

    Pulling the six into this file would make every consumer pay for all of
    them -- `hatteras_site_config` alone reads several CSVs at import -- and
    would turn the one real cycle in here (site_config imports topo_version;
    figure_style lazily imports site_config) into an import-order problem.

NOT FOR ANOTHER SITE
    `cascade_pipeline/` ships no site content by contract: a different study
    site writes its own sibling of this package and never touches it. One leak
    already exists (`cascade_pipeline/hindcast.py` imports hat_topo_version);
    folding site content into the package would make it permanently Hatteras-
    only.
"""
