# input_prep - building what the model eats

One folder per stage, numbered to match `data/hatteras_init/`. A stage's code
is here; everything it reads or writes is in the data tree under the same
number.

```
0-elevation/          the DEMs, from survey to 10 m domain rasters
1-barrier3d-domains/  the extraction: interiors and dune arrays per domain
2-brie-offset/        where each domain sits cross-shore at model year zero
3-env-forcings/       sea level, storms, waves
4-mgmt-forcings/      NC-12 and nourishment
5-scr/                the observed shoreline record and its rate fits
6-scr-smooth/         smoothing the observed rates
7-source-sink/        the background-erosion calibration
8-overwash-analysis/  modelled overwash against the record
```

Note the plural: the code folder is `4-mgmt-forcings`, the data folder is
`4-mgmt-forcing`. That is a spelling accident, not a distinction.

## Before running any of these

Most take the period as an argument and derive every path from it. Prefer that
to editing a literal: a folder name and a file name that are typed separately
are a pair that can disagree, and have.
