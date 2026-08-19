"""Rendering modules: shoreline_gif (animation) and rate_comparison (static figures).

Both import the shared foundation (domains, run_info, annotations,
coastsat_loess) but never import each other -- keep it that way. If a
helper is needed by both, it belongs one level up, not duplicated here.
"""
