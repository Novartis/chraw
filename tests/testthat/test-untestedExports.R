## Coverage for the exported functions that had no test at all.

test_that("the metadata accessors return the recorded values", {

    expect_identical( referenceGenome(ce_examples), "mm10" )
    expect_identical( referenceGenome(ce_atac), "hg38" )
    expect_identical( referenceGenome(ce_rat), "rn6" )
    expect_identical( pipeline(ce_examples), "public" )

    expect_error( referenceGenome("not a ChrawExperiment") )
})

test_that("bam files can be indexed and fragment lengths computed", {

    ce <- freshExample("ce_chipseq")
    indexBam( ce )
    bams <- colData(ce)$bamFile
    expect_true( all(file.exists(paste0(bams, ".bai"))) )

    fld <- computeFragLengthDist(
        bams[1], param = csaw::readParam(pe = "both", restrict = "chr6") )
    expect_s3_class( fld, "data.frame" )
    expect_named( fld, c("bamReads", "fraglen", "freq") )
    expect_true( all(fld$freq >= 0 & fld$freq <= 1) )
    expect_equal( sum(fld$freq), 1, tolerance = 1e-8 )

    ce <- addFragmentLengthDist(
        ce, param = csaw::readParam(pe = "both", restrict = "chr6") )
    expect_true( "FragLengthDist" %in% names(experiments(ce)) )
    mat <- experiments(ce)[["FragLengthDist"]]
    expect_identical( colnames(mat), rownames(colData(ce)) )
    expect_true( all(mat >= 0) )

    p <- plotFragmentLengthDist( ce )
    expect_s3_class( p, "ggplot" )
    expect_identical( p$labels$y, "Fraction of fragments" )
})

test_that("narrow peaks can be imported from a single file", {

    peakFile <- colData(ce_examples)$peakFile[1]
    peaks <- importNarrowPeaksFromFile( peakFile )
    expect_s4_class( peaks, "GRanges" )
    expect_true( length(peaks) > 0 )
    expect_true( all(c("name", "score", "signalValue") %in%
                     colnames(mcols(peaks))) )
    expect_error( importNarrowPeaksFromFile("does-not-exist.narrowPeak.gz") )
})

test_that("the QC plots that read the json report are labelled", {

    ## 'frac_mito' and 'frac_reads_in_annot' are only reported by the ATAC-seq
    ## pipeline, so these two plots are exercised on the ATAC example.
    p <- plotChrMStats( ce_atac )
    expect_s3_class( p, "ggplot" )
    expect_identical( p$labels$fill, "Chromosome" )
    expect_identical( p$labels$x, "Sample (replicate)" )

    p <- plotFracReadsInAnnot( ce_atac )
    expect_s3_class( p, "ggplot" )
    expect_identical( p$labels$fill, "Genomic annotation" )
    ## The legend keys must be the annotation names, not column indices.
    expect_true( any(grepl("Promoter|Enhancer|DNase", p$data$anno)) )

    ## A ChIP-seq report has neither section, and must say so plainly.
    expect_error( plotChrMStats( ce_examples ), "frac_mito" )
    expect_error( plotFracReadsInAnnot( ce_examples ), "frac_reads_in_annot" )
    expect_error( plotFripStats( ce_examples ), "FRIP statistics" )

    expect_error( plotChrMStats("not a ChrawExperiment"),
                 "must be a ChrawExperiment" )
})

test_that("a count experiment can be turned into a tidy data frame", {

    tidyDF <- tidyChrawExperiment( ce_examples, experimentName = "Peaks" )
    expect_s3_class( tidyDF, "data.frame" )
    expect_true( all(c("sampleName", "value", "RowID") %in% colnames(tidyDF)) )

    se <- experiments(ce_examples)[["Peaks"]]
    expect_identical( nrow(tidyDF), nrow(se) * ncol(se) )
    expect_setequal( unique(tidyDF$sampleName), colnames(se) )

    tidyCpm <- tidyChrawExperiment( ce_examples, experimentName = "Peaks",
                                   assayName = "cpms" )
    expect_false( identical(tidyDF$value, tidyCpm$value) )
})

test_that("an RNA-seq experiment can be added and peaks assigned to genes", {

    skip_if_not_installed("org.Mm.eg.db")

    ## A small RangedSummarizedExperiment keyed by ENSEMBL gene identifiers.
    genesGR <- GRanges(
        rep("chr6", 4),
        IRanges(start = c(67100000, 67120000, 67140000, 67160000),
                width = 1000))
    names(genesGR) <- c("ENSMUSG00000000001", "ENSMUSG00000000003",
                        "ENSMUSG00000000028", "ENSMUSG00000000037")
    counts <- matrix(seq_len(4 * ncol(ce_examples[["Peaks"]])),
                     nrow = 4,
                     dimnames = list(names(genesGR),
                                     colnames(ce_examples[["Peaks"]])))
    rnaSE <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = counts), rowRanges = genesGR)

    expect_error(
        addRNASeqExperiment( ce_examples, rnaSE, name = "RNA",
                            identifierType = "NOT_AN_ID" ),
        "ENTREZID" )
    expect_error(
        addRNASeqExperiment( "not a ChrawExperiment", rnaSE, name = "RNA" ),
        "must be a ChrawExperiment" )

    ce <- addRNASeqExperiment( ce_examples, rnaSE, name = "RNA" )
    expect_s4_class( ce, "ChrawExperiment" )
    expect_true( "RNA" %in% names(experiments(ce)) )
    expect_true( all(c("ENSEMBL", "SYMBOL", "ENTREZID") %in%
                     colnames(rowData(experiments(ce)[["RNA"]]))) )
})

test_that("the shipped example data carries a significant differential region", {

    ## Guards against a silent regression of the example data: if a future
    ## reduction removes the reads that carry the signal, this fails loudly.
    se <- experiments(ce_examples)[["Bins"]]
    rd <- rowData(se)

    expect_true( "Bins" %in% names(experiments(ce_examples)) )
    expect_true( all(c("agonist_3h_padj", "agonist_27h_padj") %in%
                     colnames(rd)) )

    tested <- sum( !is.na(rd$agonist_27h_padj) )
    expect_gt( tested, 10 )

    hits <- sum( rd$agonist_27h_padj < 0.1, na.rm = TRUE ) +
            sum( rd$agonist_3h_padj < 0.1, na.rm = TRUE )
    expect_gt( hits, 0 )
    expect_lt( min(rd$agonist_27h_padj, na.rm = TRUE), 0.05 )
})

test_that("regions without a nearby promoter are annotated instead of failing", {

    skip_if_not_installed("TxDb.Mmusculus.UCSC.mm10.knownGene")
    skip_if_not_installed("org.Mm.eg.db")
    skip_if_not_installed("BSgenome.Mmusculus.UCSC.mm10")

    ## This window is a gene desert: no region in it overlaps a promoter, which
    ## used to abort annotateExperimentRegions().
    desert <- unlist(tile(
        GRanges("chr6", IRanges(67090000, 67170000)), width = 2000))
    names(desert) <- sprintf("bin%0.5d", seq_along(desert))

    ce <- addCountExperiment( ce_examples, regions = desert, name = "Desert",
                             BPPARAM = BiocParallel::SerialParam() )

    expect_warning(
        ce <- annotateExperimentRegions( ce, "Desert", download = FALSE ),
        "overlaps a promoter" )

    rd <- rowData(experiments(ce)[["Desert"]])
    expect_true( all(is.na(rd$active_geneIds)) )
    expect_true( all(is.na(rd$active_geneSymb)) )
    ## The annotations that do not depend on active genes are still there.
    expect_false( any(is.na(rd$annotation)) )
    expect_true( "simple_annotation" %in% colnames(rd) )
})
