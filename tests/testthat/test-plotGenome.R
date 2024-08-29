test_that("Plotting genome functions work",{
  data("ce_examples")

  ce_examples <- rewrite_paths(ce_examples)

  includeSamples <- rownames(colData(ce_examples))[1:3]
  plottingRegion <- GRanges("chr6:67090000-67170000")

  expect_error(importAndAverage( object = "nonsense", includeSamples = includeSamples, plottingRegion = plottingRegion ),"must be a ChrawExperiment")
  expect_error(importAndAverage( object = ce_examples, includeSamples = "nonsense", plottingRegion = plottingRegion ),"not found in the ChrawExperiment object.")
  expect_error(importAndAverage( object = ce_examples, includeSamples = includeSamples, plottingRegion = "nonsense" ),"plottingRegion must be")
  covData <- importAndAverage( object = ce_examples, includeSamples = includeSamples, plottingRegion = plottingRegion )
  expect_type(covData,"list")

  expect_error(plotCovFromDF(covData="nonsense", replicate_group = "H3K27ac_NR1I3_CTRL_3h", plottingRegion = plottingRegion),"must be a dataframe")
  expect_error(plotCovFromDF(replicate_group = "H3K27ac_NR1I3_CTRL_3h", plottingRegion = plottingRegion),"is missing")
  expect_error(plotCovFromDF(covData=covData, replicate_group = "nonsense", plottingRegion = plottingRegion),"must be a string containing a valid")
  expect_error(plotCovFromDF(covData=covData, plottingRegion = plottingRegion),"is missing")
  expect_error(plotCovFromDF(covData=covData, replicate_group = "H3K27ac_NR1I3_CTRL_3h", plottingRegion = "nonsense"),"must be a GRanges object")

  plot <- plotCovFromDF(covData=covData, replicate_group = "H3K27ac_NR1I3_CTRL_3h", plottingRegion = plottingRegion)
  expect_s3_class(plot,"gg")
})
