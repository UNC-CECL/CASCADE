# figure_making/model — how the models work

Figures that explain Barrier3D, BRIE, the CASCADE coupling and the management
modules, drawn from a finished Hatteras run rather than a cartoon. Output goes
to `output/figures/4-model-mechanics/<model>/`.

```
model_mechanics_figures.py   Barrier3D storm year and budget, BRIE diffusion and
                             asymmetry, the coupling split and loop, management
overwash_routing_figures.py  where water and sand go in a storm; routing checks
storm_replay.py              module: replays storms through a saved domain's own
                             update() and captures the routing arrays
```

What not to trust: `storm_replay` reproduces a saved run exactly
(`storm_routing_check`), but its `defects=` variant is an in-memory patch of
the upstream Barrier3D code, not a run anyone made. Its constants (a unit and
the source-text patches) are deliberately not a CONFIG block: they are not
settings.

## The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

### model_mechanics_figures.py

How the models work, drawn from a finished Hatteras run: Barrier3D, BRIE, the coupling, management.

From the script's original header:

```text
How the models work, drawn from a finished Hatteras run rather than from a
cartoon: the Barrier3D grid through one storm year, one domain's cross-shore
budget over a window, BRIE's alongshore diffusion, how CASCADE splits a
shoreline's change between the two, the whole island as the coupled model
holds it, the annual coupling loop with the unit handed across each join, and
what the two management modules do to the grid.

    python scripts/figure_making/model/model_mechanics_figures.py [--only NAME]

Writes to output/figures/4-model-mechanics/<model>/ (PNG at the top, PDF + CAPTIONS.md under
supporting/):

    barrier3d_storm_year.png     one domain, one storm year: the year's storms
                                 against the dune, the grid before and after,
                                 and what moved
    barrier3d_domain_budget.png  one domain over 1996-2010: the profile at the
                                 two ends, the shoreface toe / shoreline /
                                 back-barrier, and the annual fluxes
    brie_diffusion.png           BRIE's angle-dependent diffusivity and where
                                 the Hatteras shoreline sits on it
    brie_domain_order.png        BRIE alone on the 1996 offset with the
                                 domains fed south->north (as run) and
                                 reversed: the order matters through the
                                 wave asymmetry
    brie_domain_orientation.png  north-up map of the numbered domains, the
                                 wave asymmetry and net drift, beside BRIE's
                                 array (index <-> GIS) in the same orientation
    brie_asymmetry_explained.png what the asymmetry counts, why its waves are
                                 head-on to positive-θ links, what that is on
                                 Hatteras, and what it does to a cape (smoothing
                                 rate, not drift; the step does not conserve sand)
    cascade_shoreline_split.png  each domain's shoreline change split into the
                                 Barrier3D cross-shore part, the source/sink
                                 (BE) part and the BRIE alongshore part
    cascade_island_grids.png     all 90 Barrier3D grids placed on the BRIE
                                 shoreline, the island as CASCADE holds it
    cascade_coupling_loop.png    the annual loop: who runs, what is handed
                                 across, and the unit it is handed in
    management_modules.png       the roadway manager (overwash cleared, dunes
                                 rebuilt) and a nourishment spreading
                                 alongshore

THE RUNS
    Everything is read from the saved Cascade object in each run's .npz
    (output/raw_runs/matrix/...), so a figure and a run cannot disagree. The
    natural run (no road, no beach/dune manager, edgeBE, 1996-2010) carries
    the physics figures because nothing human touches its grids; the
    management figure pairs a managed run with the natural run of the same
    window. The unit contract behind the loop figure is UNITS.md.
```

Notes that were in the code:

```text
The storms as the model met them: replayed through update() from the
saved grid, which reproduces the run exactly (storm_routing_check). The
crest they are tested against is the dune AFTER the year's growth,
computed once before the first storm (storm_replay.py).
```

```text
brie.py:1310 builds A row by row from that row's own r (periodic):
A x = (1 + 2r) x - r (x[i-1] + x[i+1]) = x - r lap(x).
```

```text
State t is the grid on 1 January of y0 + t; model year t (storms of
calendar year y0 + t - 1) runs between states t - 1 and t. The roadway
series are written at index t by update t (roadway_manager.py:780).
```

<details><summary>Function notes (the original docstrings)</summary>

**`plan_grid()`**

```text
(cross-shore rows, alongshore columns) in m MHW: the two dune rows on
top of the berm, then the interior. Ocean first, as Barrier3D stores it.
```

**`brie_alone()`**

```text
BRIE's alongshore diffusion on its own (no Barrier3D, no source/sink):
the padded offset added to BRIE's straight initial shoreline, run
ORDER_YEARS annual steps. Returns the change in x_s, m, landward positive.
```

**`brie_diffusivity()`**

```text
BRIE's wave-climate diffusivity (m2/yr) at shoreline angles theta, read
from its own table exactly as the solve does (brie.py:1293, before the
clamp at zero).
```

**`split_shoreline_change()`**

```text
Invert BRIE's implicit solve year by year to recover the Barrier3D
shoreline change it was handed, so each domain's total change splits
exactly into what Barrier3D did (cross-shore) and what BRIE did
(alongshore). Only valid for a run with no management: the managers move
x_s between the two solves. All in metres, + = landward.
```

</details>

### overwash_routing_figures.py

Where the water and sand go during a storm, and whether Barrier3D's overwash routing behaves as storms grow.

From the script's original header:

```text
Where the water and the sand go during a storm, and whether Barrier3D's
overwash routing does what it should as a storm gets stronger.

    python scripts/figure_making/model/overwash_routing_figures.py [--only NAME]

Barrier3D routes each storm hour by hour over the domain grid (plus one dune
row in front and a strip of bay behind), tracking discharge and sediment flux
in every cell. None of that is saved. This script REPLAYS storms through the
model's own `Barrier3d.update()` from a saved run's grid, and reads the
routing arrays out of the running update at the point the storm finishes (a
line trace on that one frame; the model code is not copied or modified).

    storm_routing_check.png     the replay of the year's real storms against
                                the grid the run saved: if they differ, the
                                replay is not the model and nothing below
                                should be trusted
    storm_routing_ladder.png    the same grid hit by four storms of rising
                                strength, collision -> run-up through the
                                gaps -> run-up over the dune -> inundation:
                                the water that crossed each cell, the
                                elevation change, and the cross-shore
                                deposition profile
    storm_routing_hours.png     one run-up storm hour by hour: where the
                                water is and what the bed has done so far

Writes to output/figures/4-model-mechanics/<model>/, ocean at the RIGHT in every plan panel.
The domain and year are the storm-year example of model_mechanics_figures
(GIS 6, the 2006 storms, natural 1996-2010 run).
```

Notes that were in the code:

```text
Two things happen to the grid after the captured line, and the saved grid
has both: update() drops trailing all-bay rows, and the year's shoreline
change drops (retreat) or adds (progradation) rows at the ocean side.
Compared unaligned, a one-cell retreat read as an 8.6 m mismatch.
```

<details><summary>Function notes (the original docstrings)</summary>

**`ladder()`**

```text
Four storms against this grid's dune: below its lowest crest, between
its lowest and mean crest, over its mean crest, and with Rlow above the
gaps (inundation). Levels are set from the grid, so the ladder means the
same thing on any domain.
```

**`extent()`**

```text
Plan extent with the dune row at x = 0 and landward positive; the axis
is then inverted so the ocean sits at the right.
```

**`fig_storm_routing_response()`**

```text
Overwash against storm strength on one grid, as the model is (all three
overwash fixes, hatteras/adopted) and with the defects put back in memory
(storm_replay.DEFECTS, the upstream code).
```

</details>

### storm_replay.py

Replay storms through a saved Barrier3D domain's own update() and read the routing arrays it discards.

From the script's original header:

```text
Replay storms through a saved Barrier3D domain's own `update()` and read the
routing arrays out of it: water discharge, sediment flux in and out, and the
elevation of every cell for every routing step, which the model computes and
throws away. Used by model_mechanics_figures.py and overwash_routing_figures.py.

The model code is not copied or modified. A line trace on the one update()
frame snapshots its locals at the first statement after each storm's routing
loop, and stops the update once the last storm is captured.

Two facts about Barrier3D that the replay depends on, both checked against a
saved run to 0.0 m (overwash_routing_figures.storm_routing_check):
  - SeaLevel() lowers DuneDomain[t - 1] IN PLACE at the start of update t, so
    a saved object's dune slice t - 1 has already lost year t's sea-level rise
    and must have it put back. The interior is lowered into a new array, so
    DomainTS[t - 1] is the true starting grid.
  - The dune crest that decides which cells overwash (DuneDomainCrest,
    barrier3d.py:1358) is computed ONCE per year, after dune growth and
    before the first storm. Every storm that year is tested against it and
    routed over it, although each storm also lowers the dune it erodes.
```

Notes that were in the code:

```text
THREE DEFECTS in how Barrier3D starts overwash, found with this replay on
2026-09-27 and FIXED in ../Barrier3D on 2026-09-28 (990c3bd, 015f11e,
e929e65; merged into hatteras/adopted, see HATTERAS_FIXES.md there). The
model now carries all three, so `fixes=` is a no-op on it and `defects=` is
the variant that matters: it puts a defect BACK into an in-memory copy of the
class (the upstream UNC-CECL code), so the fix's effect can still be measured.
"gaps"     DuneGaps() (2020) drops the last overtopped cell of the last
gap, and returns nothing at all when one cell is overtopped.
"slice"    update() sets gap water with Discharge[:, 0, start:stop], but
stop is inclusive: every gap loses its last cell of water, and
a one-cell gap gets none (2020).
"momentum" update() computes the inundation momentum constant
C = Cx * AvgSlope before the gap loop, then resets C = 0 inside
it (b11b880, 2024 Numba refactor), so inundation transport
Ki * (Q * (S + C))**mm always runs with C = 0.
```

```text
SeaLevel() lowers DuneDomain[t - 1] IN PLACE at the start of update t
(barrier3d.py:19), so the saved object's slice t - 1 has already had
year t's sea-level rise taken off. Put it back, or the replay lowers the
dune twice. The interior is lowered into a new array, so DomainTS[t - 1]
is the true starting grid. (A 4 mm slip here moved individual cells by
up to 0.5 m while the storm's total overwash changed by 1%: the routing
is cell-scale sensitive, the volumes are not.)
```

<details><summary>Function notes (the original docstrings)</summary>

**`_dune_gaps_upstream()`**

```text
DuneGaps as upstream Barrier3D has it (UNC-CECL master): drops the last
overtopped cell of the last gap, and a lone overtopped cell entirely.
```

**`model_class()`**

```text
(class, update code object, capture line). With neither, Barrier3d
itself; otherwise a subclass whose update() is compiled from the model's
own source with the named text substitutions. A fix already in the model
is skipped; a defect is put back from its upstream text.
```

**`replay()`**

```text
Run `storms` (rows of Rhigh m MHW, Rlow m MHW, period s, duration h)
through model year `t` of a saved Barrier3D domain, starting from the
grid the run saved entering that year. Returns one dict per storm.
```

</details>
