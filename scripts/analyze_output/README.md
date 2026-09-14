# analyze_output - questions that span more than one run

Nothing here runs the model. Each script reads finished runs and writes to
`output/comparisons/`.

```
compare_runs/         one run against another, or against the observed rates
overwash/             modelled overwash against the observed record
smoothing_vs_cascade/ what the LOESS smoothing does to the comparison
```

The run tree is addressed through `cascade_pipeline.run_registry`, never by
building a path: a run's NAME describes its scenario and its PATH describes its
forcing, and the two are easy to get wrong by hand.
