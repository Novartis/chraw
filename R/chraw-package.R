#' @title chraw: analysis of chromatin and multi-omics datasets
#'
#' @description
#' `chraw` analyses chromatin and multi-omic experiments. It extends the
#' [MultiAssayExperiment::MultiAssayExperiment] class into a
#' [ChrawExperiment-class] container that additionally records the reference
#' genome and the processing pipeline of the experiment, and it collapses each
#' major analysis step (quality control, read counting, differential testing,
#' peak annotation, enrichment and visualisation) into a single function call.
#'
#' The typical workflow is:
#'
#' 1. Build the container with [ChrawExperimentFromSampleFile()] or
#'    [ChrawExperimentFromDataFrame()].
#' 2. Inspect the quality of the experiments with [plotMappingStats()],
#'    [plotPCRDuplicateStats()], [plotChrMStats()], [plotFripStats()],
#'    [plotFracReadsInAnnot()] and [plotFragmentLengthDist()].
#' 3. Define the regions of interest with [importNarrowPeaks()] or
#'    [getGenomicBins()] and count reads over them with
#'    [addCountExperiment()].
#' 4. Test for differential signal with [testForDiffSignal()] and retrieve the
#'    results with [pullDiffResults()].
#' 5. Annotate the regions with [annotateExperimentRegions()] and test
#'    annotation enrichments with [enrichAnno()].
#'
#' See `vignette("chraw", package = "chraw")` for a worked example on both
#' ATAC-seq and ChIP-seq data.
#'
#' @keywords internal
#'
#' @importFrom GenomeInfoDb keepStandardChromosomes dropSeqlevels seqlengths
#'   seqlevelsStyle seqlevelsStyle<-
#' @importFrom AnnotationDbi columns
#' @importFrom GenomicFeatures asGFF genes
#' @importFrom rtracklayer export
#' @importFrom SummarizedExperiment SummarizedExperiment assay assays
#'   assayNames rowData rowData<- rowRanges rowRanges<-
#' @importFrom IRanges CharacterList IRanges NumericList
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom GenomeInfoDb seqlevels<- seqlevelsInUse
#' @importFrom DESeq2 estimateSizeFactors sizeFactors<-
#' @importFrom methods as callNextMethod is new slot validObject
## complete.cases is deliberately not imported from stats: MultiAssayExperiment
## already exports a generic of that name and dispatches to stats for
## data frames.
#' @importFrom stats aggregate ave cor fisher.test median
#'   na.omit prcomp quantile reorder reshape sd
#' @importFrom utils head read.csv read.delim read.table tail
"_PACKAGE"

## Column names that are referenced inside `ggplot2` aesthetics and `subset()`
## style calls are not visible to `R CMD check`; declare them here.
utils::globalVariables(c(
    "anno", "annoName", "baseMean", "chromHMMDict", "cond", "conf1", "conf2",
    "contrastName", "cov", "end", "frac_mito_reads", "fragmentSize",
    "frequency", "freq", "frip", "labelHeight", "log2FoldChange",
    "mapped_reads", "non_mito_reads", "odds", "padj", "pct_duplicate_reads",
    "pct_mapped_reads", "pvalue", "read_number", "relPos", "rep",
    "replicate_group", "rpm", "samp", "sampleGroup", "sampleID",
    "significant", "start", "time", "type"
))
