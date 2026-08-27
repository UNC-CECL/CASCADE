"""cascade_pipeline: post-run analysis and figures for CASCADE hindcasts.

Geometry, shoreline extraction, CoastSat-style transect/LOESS processing,
and figure/GIF rendering for a completed CASCADE run. Consumes a run's
output; it doesn't drive the simulation itself (that stays in your own run
script). Ships with no site content -- see e.g. hatteras_site_config.py for
how one study site (Hatteras Island) supplies its own domain geometry and
place names on top of this package.

Import from submodules explicitly rather than from the package root, e.g.:

    from cascade_pipeline.domains import DomainGeometry
    from cascade_pipeline.run_info import RunInfo
    from cascade_pipeline.shoreline import build_shoreline_matrix, compute_change_rate, compute_lrr
    from cascade_pipeline.coastsat_loess import CoastSatDataset, LoessConfig, build_coastsat_series
    from cascade_pipeline.plotting.shoreline_gif import GifConfig, make_all_shoreline_gifs
    from cascade_pipeline.plotting.rate_comparison import plot_rate_comparison, plot_annotated_rate_comparison
"""
