
test_that("Test getGenomicBins dispatch", {
    data("ce_chipseq")

    ce_examples <- rewrite_paths(ce_examples)

    expect_error( selectBSgenome( "other" ), "Input needs to be")
    expect_s4_class( selectBSgenome( ce_chipseq ), "BSgenome" )
    expect_s4_class( selectOrgDb( ce_chipseq ), "OrgDb" )
    expect_s4_class( selectTxDb( ce_chipseq ), "TxDb" )
    expect_message( rs <- selectChromHMMData( ce_chipseq ), "As no ChromHMM file was specified")
    expect_type( rs, "character" )

    expect_error( getGenomicBins(ce_chipseq), "is missing, with no default" )
    expect_error( getGenomicBins(ce_chipseq, binSize="a"), "must be a numeric value" )
    expect_equal( getGenomicBins(ce_chipseq, binSize=100000), getGenomicBins(selectBSgenome(ce_chipseq), binSize=100000) )

    ce_chipseq@referenceGenome <- "rn6"
    expect_warning( x <- selectChromHMMData( ce_chipseq ), "No default ChromHMM data was found" )
    expect_null( x )
})


test_that("Test c function dispatch", {
    data(ce_rat)
    data(ce_chipseq)
    data(ce_atac)

    ce_rat <- rewrite_paths(ce_rat)
    ce_chipseq <- rewrite_paths(ce_chipseq)
    ce_atac <- rewrite_paths(ce_atac)


    ce_rat@referenceGenome <- "hg38"
    ce_atac@referenceGenome <- "hg38"
    ce_chipseq@referenceGenome <- "hg38"

    ce1 <- ce_atac
    ce2 <- ce_chipseq
    ce3 <- ce_rat

    #experiments(ce1)[2:4] <- NULL
    #experiments(ce2)[2:4] <- NULL
    experiments(ce3)[2] <- NULL

    bSampleNames <- sort(c(rownames(colData(ce1)), rownames(colData(ce2))))
    ceMerged <- c(ce1, ce2)
    expect_s4_class( ceMerged, "ChrawExperiment")
    aSampleNames <- sort(rownames(colData(ceMerged)))
    expect_equal( aSampleNames, bSampleNames )

    ceMerged <- c( ce3, ce1, ce2)
    bNames <- sort(c(rownames(colData(ce3)), rownames(colData(ce1)), rownames(colData(ce2))))
    aNames <- sort(rownames(colData(ceMerged)))
    expect_s4_class( ceMerged, "ChrawExperiment" )
    expect_equal( aNames, bNames )

    ce2@referenceGenome <- "mm10"
    expect_error(c(ce1, ce2), "different species is not supported")

    ce2@referenceGenome <- "hg38"
    ce2@pipeline <- "chromatinalign"
    expect_warning(c(ce1, ce2), "from different pipelines")

} )
