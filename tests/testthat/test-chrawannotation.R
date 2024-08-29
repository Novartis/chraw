
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

  enrich.res <- enrichAnno( ce_examples, 'Peaks', 'test', 'chromHMM_annotation_simple')
  #expect_true(is.data.frame(enrich.res))
  #expect_true(sum(colnames(enrich.res) %in% c('annoName', 'n_anno_fg_peaks',
  #                                            'n_anno_bg_peaks', 'conf1', 'conf2',
  #                                            'odds', 'pvalue')) == 7)

  expect_error(plotEnrichResults(enrich.res[,-1]), 'Input data.frame is missing')
  #expect_s3_class( plotEnrichResults(enrich.res), "gg" )
})
