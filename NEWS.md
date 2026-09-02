# chraw 0.99.0

* First submission to Bioconductor.
* `ChrawExperiment` gains the accessors `referenceGenome()` and `pipeline()`;
  code should no longer reach into the slots directly.
* `annotateExperimentRegions()` now defaults to `download = FALSE`. ChromHMM
  annotation files that are downloaded are cached with `BiocFileCache` instead
  of being re-fetched on every call.
* All plotting functions now label every axis, legend and colour scale with the
  quantity being shown and its unit, and share a common theme.
* `plotPCA()` no longer duplicates the word "variance" in its axis titles.
* `plotChrMStats()` now reports fragments on the same scale as
  `plotMappingStats()`, and `plotFripStats()` reports the same unit as
  `plotFracReadsInAnnot()`.
* `plotVProfile()` labels its colour scale according to the `rowScaled`
  argument, and its `sampleLabelColumn = NULL` default now works.
* `plotDiffScatter()` names its colour legend after the `colorByFDR` argument
  and no longer emits `NA` in axis titles.
* The example data shipped in `inst/files` and `data/` was reduced to the
  genomic windows exercised by the vignette and the unit tests.
