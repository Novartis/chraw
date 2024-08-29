
test_that("QC plotting functions work", {
    data("ce_chipseq")
    ce_chipseq <- rewrite_paths(ce_chipseq)

    experiments(ce_chipseq)[2] <- NULL
    plot_ms <- plotMappingStats( ce_chipseq )
    expect_s3_class( plot_ms, "gg" )
    expect_error( plotMappingStats( "wronginput" ), "must be a ChrawExperiment object")
    plot_pcr <- plotPCRDuplicateStats( ce_chipseq )
    expect_s3_class( plot_pcr, "gg" )
    expect_error( plotPCRDuplicateStats("wronginput"), "must be a ChrawExperiment" )
    expect_error( plotFragmentLengthDist("wronginput"), "must be a ChrawExperiment" )
    expect_error( plotFragmentLengthDist(ce_chipseq), "Fragment length distribution" )
    data("ce_examples")
    ce_examples <- rewrite_paths(ce_examples)
    colnames(colData(ce_examples))[colnames(colData(ce_examples)) == "sampleName"] <- "error"
    expect_error( plotFragmentLengthDist(ce_examples), "should be the name of the samples column." )
} )

test_that( "Exploratory plots work as expected", {
    data(ce_examples)

    ce_examples <- rewrite_paths(ce_examples)

    expect_error(plotPCA( ce_examples, experimentName="nonsense" ), "was not found in the")
    expect_error(plotPCA( ce_examples, "Peaks", colourGrouping="nonsense" ), "not found in the ChrawExperiment")
    expect_error(plotPCA( ce_examples, "Peaks", prinComps=c(1, Inf) ), "out of range")
    xx <- plotPCA( ce_examples, "Peaks", colourGrouping="condition" , transformation=asinh)
    expect_s3_class(xx, "gg")
    sampleSubset <- rownames(colData(ce_examples))[1:2]
    xx <- plotPCA( ce_examples, "Peaks", colourGrouping="condition", includeSamples=sampleSubset ,transformation=asinh)
    expect_true(all(xx$data$condition %in% colData(ce_examples)[sampleSubset,"condition"]))
    expect_error(plotPCA( ce_examples, "Peaks", colourGrouping="condition", includeSamples=c(sampleSubset, "sheet1") ),
                 "not found")

    xx <- plotSampleCorrelations( ce_examples, "Peaks", annotationColumns="condition" )
    expect_s4_class(xx, "Heatmap" )

    expect_error(plotSampleCorrelations( ce_examples, experimentName="nonsense" ), "was not found in the")
    expect_error(plotSampleCorrelations( ce_examples, "Peaks", annotationColumns="nonsense" ), "not found in the ChrawExperiment")
    expect_s4_class( plotSampleCorrelations( ce_examples, "Peaks"), "Heatmap" )

    xx2 <- plotSampleCorrelations(
        ce_examples, "Peaks", annotationColumns="condition",
        annotationColors=list(condition=c(WIZ="red", `WIZ_GenScript_11`="blue", `WIZ_Genscript_15`="green")),
        show_row_names = FALSE, show_column_names = FALSE )
    expect_s4_class(xx2, "Heatmap")

    expect_error(plotSampleCorrelations(
        ce_examples,
        experimentName="Peaks", annotationColumns="condition",
        includeSamples=c("H3K27ac_NR1I3_CTRL_3h_rep1",
                         "H3K27ac_NR1I3_CTRL_27h_rep2", "nonsense")),
        "not found in the ChrawExperiment object" )

    colData(ce_examples)$sampleAlias <- paste0("sample", seq_len(nrow(colData(ce_examples))))

    xx3 <- plotSampleCorrelations(
        ce_examples,
        experimentName="Peaks", annotationColumns="condition",
        includeSamples=c("H3K27ac_NR1I3_CTRL_3h_rep1",
                         "H3K27ac_NR1I3_CTRL_3h_rep2"),
        sampleLabelColumn="sampleAlias" )

    expect_s4_class(xx3, "Heatmap")

} )

test_that("Functions to get metaprofiles work well", {
    data("ce_examples")

    ce_examples <- rewrite_paths(ce_examples)

    wizExperiments <- rownames(colData(ce_examples))[1:2]
    peaks <- importNarrowPeaks(
        ce_examples, includeSamples=wizExperiments, merge=TRUE )
    names(peaks) <- sprintf("peak%0.6d", seq_len(length(peaks)))
    rgs <- resize( peaks, width=3000, fix="center" )
    rgs <- rgs[seqnames(rgs) == 'chr6' & start(rgs) >= 67090000 & end(rgs) <= 67170000]
    expect_error(calculateMetaProfiles( 3, regions=rgs ), "is not TRUE")
    expect_error(calculateMetaProfiles( ce_examples, regions=rgs, inputColumn="nonsense" ),
                 "Invalid `inputColumn=` parameter")
    expect_error(calculateMetaProfiles( ce_examples, regions=rgs, includeSamples="nonsense" ),
                 "not found in the ChrawExperiment object" )
    expect_error(calculateMetaProfiles( ce_examples, regions=rgs, outputFormat="nonsense" ),
                 "The parameter 'outputFormat' must be either" )
    rs1 <- calculateMetaProfiles( ce_examples, regions=rgs, outputFormat="ScoreMatrixList" )
    expect_s4_class( rs1, "ScoreMatrixList")
    expect_identical(nrow(colData(ce_examples)), length(rs1))
    rs2 <- calculateMetaProfiles( ce_examples, regions=rgs, includeSamples=wizExperiments, outputFormat="list" )
    rs3 <- reformatMatrixScores( rs1, ce_examples, outputFormat="list" )
    expect_identical( rs3[wizExperiments], rs2[wizExperiments] )

    rs4 <- reformatMatrixScores( rs1, ce_examples, outputFormat="long_df" )
    expect_identical(
        rs3[[wizExperiments[1]]]["peak055250",],
        rs4$score[rs4$sampleName == wizExperiments[1] & rs4$regionID == "peak055250"] )

    expect_error(
        reformatMatrixScores( rs1, object=ce_examples, outputFormat="long_df", annotationColumns="nonsense"),
        "not found in the ChrawExperiment object" )

    expect_error(
        plotMetaProfileHeatmaps( ce_examples, rgs, includeSamples=wizExperiments, sampleLabelColumn="nonsense" ),
        "not found in the ChrawExperiment object")

    colData(ce_examples)$sampleAlias <- colData(ce_examples)$sampleName
    colData(ce_examples)$sampleAlias <- paste0("sample", seq_len(nrow(colData(ce_examples))))

    ht <- plotMetaProfileHeatmaps( ce_examples, rgs, includeSamples=wizExperiments,
                                  sampleLabelColumn="sampleAlias", bin.num=100 )
    expect_s4_class(ht, "HeatmapList")

    vprof <- plotVProfile( ce_examples, rgs, sampleLabelColumn="sampleAlias" )
    expect_s3_class(vprof, "ggplot")
    expect_error(
        plotVProfile( ce_examples, rgs, sampleLabelColumn="sampleAlias2" ),
        "not found in the ChrawExperiment")

    expect_error(
        plotVProfile( ce_examples, rgs, sampleLabelColumn="sampleAlias2", includeSamples="camara"),
        "not found in the ChrawExperiment")

    xx <- plotVProfile( ce_examples, rgs, sampleLabelColumn="sampleAlias", returnMatrix=TRUE )
    expect_s3_class(xx, "data.frame")

    head( xx )

    xx$rpm <- 10^6 * xx$n/xx$totReads

    # xx2 <- xx %>%
    #     dplyr::group_by( samp, fragmentSize ) %>%
    #     dplyr::mutate( rpm=(rpm-mean(rpm, na.rm=TRUE)))

    xx2 <- transform(xx, rpm = ave(rpm, samp, fragmentSize, FUN = function(x) x - mean(x, na.rm = TRUE)))

    xx <- transform(xx, rpm = ave(rpm, samp, fragmentSize, FUN = function(x) (x - mean(x, na.rm = TRUE))/sd(x, na.rm=TRUE)))
    # xx <- xx %>%
    #     dplyr::group_by( samp, fragmentSize ) %>%
    #     dplyr::mutate( rpm=(rpm-mean(rpm, na.rm=TRUE))/sd(rpm, na.rm=TRUE) )


        head(xx2[which(is.na(xx2$rpm)),])

    vprof <- plotVProfile( ce_examples, rgs, sampleLabelColumn="sampleName", rowScaled=TRUE )
    expect_s3_class(vprof, "ggplot")

} )
