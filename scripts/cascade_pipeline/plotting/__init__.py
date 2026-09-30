"""Rendering modules: shoreline_gif (animation) and rate_comparison (static figures).

Both import the shared foundation (domains, run_info, annotations,
coastsat_lowess) but never import each other -- keep it that way. If a
helper is needed by both, it belongs one level up, not duplicated here.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""
