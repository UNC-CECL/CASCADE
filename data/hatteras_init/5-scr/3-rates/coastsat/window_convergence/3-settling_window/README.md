# 3-settling_window — how many years does each place need?

**How many years of record does each place need before its shoreline change
rate stays within X m/yr of the 1996–2024 rate?** One figure answers it:

![years needed](years_needed_alongshore.png)

`years_needed_alongshore.png` has one line of numbers per transect (906) along
the island, south to north. (a) is windows starting 1996, with the end moving
later, and (b) is windows ending 2024, with the start moving earlier. The shaded
bands are three thresholds, ±1.0, ±0.5 and ±0.25 m/yr: the top of each band is
the years that threshold needs.

- **(a) forward**: 15 yr is enough at: 36% (±1.0 m/yr) · 10% (±0.5 m/yr) · 0% (±0.25 m/yr) of transects
- **(b) backward**: 15 yr is enough at: 46% (±1.0 m/yr) · 19% (±0.5 m/yr) · 5% (±0.25 m/yr) of transects

"Stays within" means the window's rate is inside the threshold and every
longer window's rate is too. The first time a rate passes through the band
doesn't count, because most transects swing through it and back out.
29 years means only the full record matches.

```
forward_from_1996/     windows 1996–2000 … 1996–2024
backward_from_2024/    windows 2020–2024 … 1996–2024
    a-eight_sites/     shoreline position through time at 8 transects, with
                       the fits drawn on it (supporting figure) + tables
    b-every_transect/  per-transect tables (the figure's source)
    c-domain_means/    the same, averaged to the 90 model domains (tables)
```

Redraw without refitting:
`python scripts/input_prep/5-scr/3-rates/coastsat/window_convergence/coastsat_window_convergence.py --figure-only`

**Retired 2026-09-29** (Hannah: too many figures, an abstract "rate minus
reference" axis, and tolerance jargon). The figures that went are
`window_convergence_*` (eight sites and domains, the error curves),
`convergence_alongshore_*`, `tolerance_comparison_*` and
`domain_mean_vs_transect_*`, in both directions. Their tables are all kept.
PNGs are gitignored, but each one's PDF was committed in 3c6de274, under the old
layout: `window_convergence/record_<start>_<end>/<direction>/{sites,all_transects,domain_means}/supporting/`.
To see one: `git show 3c6de274:data/hatteras_init/5-scr/3-rates/coastsat/window_convergence/record_1996_2024/forward_from_1996/all_transects/supporting/convergence_alongshore_forward_from_1996.pdf > old.pdf`.
