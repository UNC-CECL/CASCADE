# 2026-10-21 groin spread check

`cascade/groin.py` gained `spread_weights` (default `(1, 0, 0)`, no spread): the updrift pad's yearly loss across
the groin face shared with the two pads beyond it. Does the default leave the groin unchanged, and does set E
`(0, 2/3, 1/3)` built in reproduce the experiment patch (`2026-10-20-spread-relaxation/spread_relaxation.py`)?

Four runs, 1996-2009, edgeBE, full management, the runner as committed at 1d36cf20: `old` (the module before the
change, strength 0.16), `new` (the new module, default, 0.16), `oldE` (the old module with the experiment patch,
set E, 0.10), `newE` (the new module with set E built in, 0.10).

| check | result |
|---|---|
| new against old, default | shoreline matrix bit-identical, groin_diagnostics.csv byte-identical |
| newE against oldE, set E | shoreline matrix bit-identical, groin_diagnostics.csv byte-identical |
| newE against the spread study's `current-wE/s0.10` run (00:07) | bit-identical |
| unit tests (tests/test_groin.py) | 17 of 17 pass (6 new: the spread in both orientations, off while the groin is, the weights' check) |
| stress rig, default against the old module | 219 runs bit-identical; 78 stopped the same way |
| stress rig, set E built in against the patch | 1,650 model years bit-identical, the spread active in 962; 6 stopped the same way |

The stress rig covers the updrift pad on the higher-index side (Buxton's); the unit tests cover both sides.
