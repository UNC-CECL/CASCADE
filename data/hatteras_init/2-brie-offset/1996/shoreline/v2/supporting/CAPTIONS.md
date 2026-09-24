# Captions — v2

Written by the figure scripts through `hat_figure_style.caption()`; the images carry no titles or footnotes, this file does.

**`Island_Shoreline_Offsets_1996_buffer_diagnostic.png`.** The padded mean shoreline offset for the 1996 start, in metres as written to Island_Shoreline_Offsets_1996_PADDED_120.csv. (a) Offset along the padded domains; the real reach GIS 1-90 in blue, the 15 buffer domains each side shaded. The buffers close BRIE's periodic domain from the last real domain back round to the first along a cubic Hermite matched to the island's end slopes (cascade_pipeline.hindcast.pad_offset_ring), the same array offset_mode 'metres' hands Cascade. (b) The shoreline angle between each padded domain and the next (atan2 of the offset step over 500 m), the last point being the wrap; dashed at +/-42 degrees, past which BRIE's alongshore diffusivity changes sign.
