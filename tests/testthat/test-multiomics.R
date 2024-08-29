test_that("Test for multiomics integration functions", {

    # Load test data
    data(ce_examples)

    ce_examples <- rewrite_paths(ce_examples)

    expect_error(assignPeaks2Peaks(object = ce_examples,
                                   experimentNames = c("Peaks", "peaksIdr2")),
        "No experiment called peaksIdr2 in" )

    expect_error(assignPeaks2Peaks(object = ce_examples,
                                   experimentNames = c("ExistingSampleFlag", "Peaks")),
        "ExistingSampleFlag in ChrawExperiment object is not a")

    expect_error(assignPeaks2Peaks(object = ce_examples,
                                   experimentNames = c("Peaks", "Peaks"),
                                   filterAnnoCols = c("test", "simple_annotation")),
                 "test is not existing in the rowData")

    expect_error(assignPeaks2Peaks(object = ce_examples,
                                   experimentNames = c("Peaks", "Peaks"),
                                   filterAnnoCols = c("simple_annotation", "simple_annotation")),
                 "Please specify a filterValue to subset the peaks for simple_annotation.")

    expect_error(assignPeaks2Peaks(object = ce_examples,
                                   experimentNames = c("Peaks", "Peaks"),
                                   filterAnnoCols = c("simple_annotation", "simple_annotation"),
                                   filterValues = c("Promoter", "test")),
                 "test is not part of the filter colunmn simple_annotation.")

    expect_error(assignPeaks2Peaks(object = ce_examples,
                                   experimentNames = c("Peaks", "Peaks"),
                                   filterAnnoCols = c("simple_annotation", "simple_annotation"),
                                   filterValues = c("Promoter", "Intron"),
                                   assignFunc = 'test'),
                 "The peak to peak assignment function is not a valid function.")

    expect_error(assignPeaks2Peaks(object = ce_examples,
                                   experimentNames = c("Peaks", "Peaks"),
                                   filterAnnoCols = c("simple_annotation", "simple_annotation"),
                                   filterValues = c("Promoter", "Intron"),
                                   assignFunc = distance),
                 "The peak to peak assignment function is not a valid function.")

    expect_s4_class(assignPeaks2Peaks(object = ce_examples,
                      experimentNames = c("Peaks", "Peaks"),
                      filterAnnoCols = c("simple_annotation", "simple_annotation"),
                      filterValues = c("Promoter", "Intron"),
                      assignFunc = nearest), class = 'DataFrame')

    expect_equal(nrow(assignPeaks2Peaks(object = ce_examples,
                      experimentNames = c("Peaks", "Peaks"),
                      filterAnnoCols = c("simple_annotation", "simple_annotation"),
                      filterValues = c("Promoter", "Intron"),
                      assignFunc = nearest)), expected = 36402)

    expect_equal(nrow(assignPeaks2Peaks(object = ce_examples,
                                        experimentNames = c("Peaks", "Peaks"),
                                        filterAnnoCols = c("simple_annotation", "simple_annotation"),
                                        filterValues = c("Distal Intergenic", "Promoter"),
                                        assignFunc = nearest)), expected = 72521)

    ce_examples <- testForDiffSignal(
      ce_examples, experimentName="Peaks",
      design=~condition,
      contrasts=list(CTRL_3h=c("condition","CTRL_3h","agonist_3h")))

    expect_error(plotDiffScatter(object = ce_examples,
                                contrasts =  c('agonist_3h', 'CTRL_3h', 'test'),
                                experimentNames = c('Peaks', 'Peaks'),
                                filterAnnoCols = c('simple_annotation', 'simple_annotation'),
                                filterValues = c('Promoter', 'Distal Intergenic'),
                                assignFunc = GenomicRanges::nearest),
                 'Too many arguments given to contrasts. Please provide only two contrasts.')

    expect_error(plotDiffScatter(object = ce_examples,
                                contrasts =  c('agonist_3h', 'CTRL_3h'),
                                experimentNames = c('Peaks', 'Peaks', 'test'),
                                filterAnnoCols = c('simple_annotation', 'simple_annotation'),
                                filterValues = c('Promoter', 'Distal Intergenic'),
                                assignFunc = GenomicRanges::nearest),
                 'Too many arguments given to experimentNames. Please provide only two experimentNames')

    expect_error(plotDiffScatter(object = ce_examples,
                                 contrasts =  c('agonist_3h', 'test'),
                                 experimentNames = c('Peaks', 'Peaks'),
                                 filterAnnoCols = c('simple_annotation', 'simple_annotation'),
                                 filterValues = c('Promoter', 'Distal Intergenic'),
                                 assignFunc = GenomicRanges::nearest),
                 'Contrast test is not available for experiment Peaks.')

    expect_true(is.ggplot(plotDiffScatter(object = ce_examples,
                                 contrasts =  c('agonist_3h', 'CTRL_3h'),
                                 experimentNames = c('Peaks', 'Peaks'),
                                 filterAnnoCols = c('simple_annotation', 'simple_annotation'),
                                 filterValues = c('Promoter', 'Distal Intergenic'),
                                 assignFunc = GenomicRanges::nearest)))

    expect_true(is.ggplot(plotDiffScatter(object = ce_examples,
                                          contrasts =  c('agonist_3h', 'CTRL_3h'),
                                          experimentNames = c('Peaks', 'Peaks'),
                                          assignFunc = GenomicRanges::findOverlaps)))
})
