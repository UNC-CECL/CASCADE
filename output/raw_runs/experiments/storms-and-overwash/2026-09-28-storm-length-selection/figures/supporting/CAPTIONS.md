# Captions — figures

Written by the figure scripts through `hat_figure_style.caption()`; the images carry no titles or footnotes, this file does.

**`stage1_scores.png`.** Storm-length candidates against the observed overwash record and CoastSat. Every candidate keeps every storm event and trims those longer than L hours to the L hours around their peak; 'full' trims nothing (a 240 h limit). Thin horizontal lines: the committed series, which drops events over 72 h. Red 1996-2010, blue 2010-2024; solid managed (full_management), dashed natural. Overwash scores use every image x domain cell the imagery assessed; a model year's overwash is shared among the storms that reached a dune gap and dated by their end (7-day grace). Edge rates as the matrix (solved on the committed series), so RMSE is not yet a fair comparison.

**`stage1_per_image.png`.** Domains with overwash in each image (grey bars, observed; only domains the image assessed) against the model's count for the same image window, managed runs, for the committed series and three candidates. Timing agreement is what the per-image correlation in stage1_scores measures.
