test_that("Functions work for other species", {

    data("ce_rat")

    inst_dir <- system.file("files/", package="chraw")
    ce_rat$bamFile <- gsub("inst/files", inst_dir, ce_rat$bamFile)
    ce_rat$bigwigFile <- gsub("inst/files", inst_dir, ce_rat$bigwigFile)
    ce_rat$peakFile <- gsub("inst/files", inst_dir, ce_rat$peakFile)
    ce_rat$qcFile <- gsub("inst/files", inst_dir, ce_rat$qcFile)

    expect_s4_class(getGenomicBins( ce_rat, binSize=10000 ), "GRanges")

    expect_s4_class( ce_rat, "ChrawExperiment" )
    expect_warning(
      ce_rat <- annotateExperimentRegions( ce_rat, "Peaks" ),
        "No default ChromHMM data was found" )
    elementTypes <- elementMetadata(mcols(rowRanges(experiments(ce_rat)[["Peaks"]])))$type
    expect_true(all(c("Chraw2 annotation", "UCSC link") %in% elementTypes))

    xx <- require(TxDb.Rnorvegicus.UCSC.rn6.refGene)
    if( xx ){
        proms <- keepStandardChromosomes(promoters(TxDb.Rnorvegicus.UCSC.rn6.refGene), pruning.mode="coarse")
        notUniqueProms <- proms
        names(notUniqueProms) <- sprintf("prom%0.4d", seq_len(length(notUniqueProms)))
        proms <- unique(proms[seqnames(proms) == "chr1"])
        proms <- head(proms, 1000)
        names(proms) <- sprintf("prom%0.4d", seq_len(length(proms)))
        expect_s4_class( ce_rat, "ChrawExperiment" )
        expect_s4_class(getGenomicBins( ce_rat, binSize=10000 ), "GRanges")
        expect_error(
            addCountExperiment( ce_rat, regions=notUniqueProms, name="perPeak2" ),
            "must be unique")
        expect_warning(
            ce_rat <- annotateExperimentRegions( ce_rat, "Peaks" ),
            "No default ChromHMM")
        elementTypes <- elementMetadata(mcols(rowRanges(experiments(ce_rat)[["Peaks"]])))$type
        expect_true(all(c("Chraw2 annotation", "UCSC link") %in% elementTypes))
    }

})
