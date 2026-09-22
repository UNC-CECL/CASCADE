Imports
Fixed dune/topo (period-independent) — plot to QC
Domain orientation: set START_YEAR, resolve period config, shoreline_offset init — plot initialization
roadway_manager: load setbacks + historical events, apply road-exists-AND-not-community rule — plot
beach_dune_manager: nourishment schedule + overwash filter — plot
hard_structures/groin: Buxton groin config — plot deterioration curve
Forcings: storms + RSLR — plot trend + storm catalog
Source/sink (BE rates) — plot base vs. calibrated
CoastSat target rates — plot LOESS windows
Plotting functions (imported from hat_pipeline.plotting, not redefined inline)
Define run_cascade_simulation()
Initialize Cascade (single config, no sweep)
Run the loop, then generate the position-over-time GIF from the saved shoreline_m

