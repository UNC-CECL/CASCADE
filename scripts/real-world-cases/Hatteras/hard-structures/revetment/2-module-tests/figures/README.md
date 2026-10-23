# figures — how the simple revetment works

    python revetment_module_schematic.py

Draws only; reads no model output.

- `revetment_how_it_works.png`: four cross-sections. No wall (dune and
  shoreline move together), wall built (dune line fixed, crest floored, gap =
  the beach), years later (gap shrinks, overwash stops), gap 0 (the repaired wall holds
  the shoreline, with sand added each year).
- `revetment_yearly_loop.png`: what the module does each model year around
  `cascade.update()`.

PDFs and captions in `supporting/`.
