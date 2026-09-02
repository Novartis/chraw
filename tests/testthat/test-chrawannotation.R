
test_that("Test validity of inputs for annotation", {
  data(ce_examples)
  ce_examples <- rewrite_paths(ce_examples)

  expect_error(
      annotateExperimentRegions( ce_examples, "Peaks2" ),
      "was not found in the ChrawExperiment" )
  expect_error(
      annotateExperimentRegions( ce_examples ),
      "missing, with no default" )
  expect_error(
      annotateExperimentRegions( ce_examples, "Peaks", activeGeneList="BRD4" ),
      "is not an Entrez ID")
  suppressWarnings(ce2 <- annotateExperimentRegions( ce_examples, "Peaks" ))
  elementTypes <- elementMetadata(mcols(rowRanges(experiments(ce2)[["Peaks"]])))$type
  expect_true(all(c("Chraw2 annotation", "UCSC link") %in% elementTypes))
} )



test_that("Test annotation enrichment function", {
  data(ce_examples)

  ce_examples <- rewrite_paths(ce_examples)

  expect_error(
    enrichAnno( ce_examples, 'Peaks', 'test', 'simple_anno', foreground = 'up',
                background = 'up'),
    "Foreground condition can't be" )

  expect_error(
    enrichAnno( ce_examples, 'Peaks2', 'test', 'simple_anno'),
    "was not found in the ChrawExperiment.")

  expect_error(
    enrichAnno( ce_examples, 'Peaks', 'test', 'simple_anno'),
    "Diff peak analysis for ")

  ce_examples <- testForDiffSignal(ce_examples, experimentName="Peaks",
                          design=~condition,
                          contrasts=list(test=c("condition", "agonist_3h", "CTRL_3h")))

  suppressWarnings( ce_examples <- annotateExperimentRegions( ce_examples, "Peaks" ) )

  ## The ChromHMM columns are only added when a ChromHMM file is supplied or
  ## downloaded, which annotateExperimentRegions() no longer does by default.
  expect_error(
    enrichAnno( ce_examples, 'Peaks', 'test', 'chromHMM_annotation_simple'),
    "Annotation column is not existing" )

  ## enrichAnno() returns NULL when nothing passes the thresholds, which is the
  ## case for the reduced example data.
  enrich.res <- enrichAnno( ce_examples, 'Peaks', 'test', 'simple_annotation')
  expect_true( is.null(enrich.res) || is.data.frame(enrich.res) )

  ## plotEnrichResults() must reject a frame that lost a required column.
  fakeRes <- data.frame(
    annoName = c("Promoter", "Intron"),
    n_anno_fg_peaks = c(5, 3), n_anno_bg_peaks = c(50, 30),
    conf1 = c(-1, -2), conf2 = c(2, 1), odds = c(0.5, -0.5),
    pvalue = c(0.01, 0.4) )
  expect_s3_class( plotEnrichResults(fakeRes), "gg" )
  expect_error(plotEnrichResults(fakeRes[,-1]), 'Input data.frame is missing')
})
