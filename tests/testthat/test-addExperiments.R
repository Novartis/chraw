
test_that("Test validity of inputs for addCountExperiment", {
    data(ce_chipseq)
    ce_chipseq <- rewrite_paths(ce_chipseq)

    bins <- getGenomicBins( ce_chipseq, binSize=100000 )
    expect_error( addCountExperiment( ce_chipseq, bins, name="countsPrueba" ), "needs to be a named GRanges object" )
    names(bins) <- sprintf("bin%0.6d", seq_len(length(bins)))
    expect_error( addCountExperiment( ce_chipseq, bins, name=2 ), "must be a character" )
    expect_error( addCountExperiment( ce_chipseq, includeSamples="notexistingsample", name="countsPrueba" ), "could not be found" )
    names(bins) <- names(bins[1])
    names(bins)[2] <- names(bins)[1]
    expect_error( addCountExperiment( ce_chipseq, bins, name="countsPrueba" ), "must be unique" )
} )

test_that( "Test validity of inputs for importNarrowPeaks", {
    data(ce_chipseq)
    ce_chipseq <- rewrite_paths(ce_chipseq)
    expect_error(importNarrowPeaks(ce_chipseq, includeSamples="notexistingsample"), "could not be found")
    expect_error(importNarrowPeaks(ce_chipseq, peakType="camara"), "Invalid 'peakType' parameter")
} )
