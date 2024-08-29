#' Annotate a set of regions (peaks)
#'
#' This function annotates the peak list of an experiment in a ChrawExperiment
#' object. Annotations are done according to their annotation with gene region
#' overlaps, annotation to the nearest TSS, annotation of genes flanking the
#' regions, annotations with gene regions considering only active genes, and
#' annotation with overlap with ChromHMM colors.
#'
#' @param object A ChrawExperiment object.
#' @param experimentName The name of the count experiment to test, normally added using either the `addCountExperiment` function.
#' @param chromHmmAnnotation Path to an annotation file (.bed format). It should correspond to the type of tissue and genome assembly used to produce the data. Annotation files can be found on \href{https://www.encodeproject.org}{ENCODE} or on the \href{https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html}{Roadmap Epigenomics Project webpage}. Annotations not using ChromHMM 15 and 18 state models may lead to incomplete simplified annotations. Please note that the annotation step will order the annotations by their name and use this order when finding overlapping peak. As such, it is advised to use ordered prefixes (eg, "1_TssA", "2_TssAFlnk") in front of the annotation names, so that they can be ordered by importance. Most chromHMM files already contain annotations and colors columns. If no file is provided, a default file will automatically be selected (NHEK cell line, hg19 and hg38 assemblies, liver cells for mm10 assembly). No default file is available for rn6.
#' @param flankDistance Flanking distance to peak used during the annotation step. Default 50Kb.
#' @param download Boolean for downloading chromHmm data. Default FALSE.
#' @param activeGeneList A list of Entrez identifiers with 'active' genes. By default (\code{NULL}), peaks are annotated as active genes if they overlap with promoter regions.
#' @return A ChrawExperiment object. Annotation columns are included in `rowData(experiments(object)[[experimentName]])`. The most relevant columns are:
#' \itemize{
#' \item{'SYMBOL', 'geneId': Symbol and gene IDs of the nearest gene that overlaps with each peaks.}
#' \item{'nearest_TSS', 'nearest_TSSIds': Symbol and gene IDs if the nearest TSS.}
#' \item{'flank_geneSymb', 'flank_geneIds': Symbol and gene IDs of all genes flanking your peak (distance controlled by 'flankDist' option). Multiple genes can be returned for the same peaks. If so, the respective symbols/IDs are ";" separated in the same cell.}
#' \item{'active_geneIds', 'active_geneSymb': Symbol and gene IDs of the active genes, defined by an overlap of the peaks with promoter regions OR by overlap with the regions passed to 'activeGeneList' option. The 'active gene' annotation is only useful for activation marks such as ATAC or TF-ChIP but NOT for methylation marks.}
#' \item{'annotation', 'simple_anno': gene-centric annotation and simplified annotations of the peaks.}
#' \item{'custom_annot', 'simple_custom_annot': annotation from the custom input that can be passed to 'annot_file' option. Usually this will be chromHMM annotation (this is the case with the default annotation file).}
#' }
#' @importFrom rtracklayer import
#' @seealso annotatePeak, annotatePeakInBatch
#' @examples
#' data(ce_examples)
#' ce_examples <- rewrite_paths(ce_examples)
#'
#' ce_examples <- annotateExperimentRegions( ce_examples, "Peaks" )
#'
#' @export
annotateExperimentRegions <- function( object, experimentName, chromHmmAnnotation = NULL,
                                      flankDistance = 50000, activeGeneList = NULL, download = TRUE ){
    validObject(object)
    if( !is(object, "ChrawExperiment") ){
        stop("The `object=` parameter must be a ChrawExperiment object")
    }
    if( !experimentName %in% names(experiments(object)) ){
        stop(sprintf("The experiment %s was not found in the ChrawExperiment object", experimentName))
    }
    ## select the appropiate metadata ##
    message("Loading annotation data...")
    bsGenome <- selectBSgenome( object )
    orgDb <- selectOrgDb( object )
    txDb <- selectTxDb( object )
    txDb <- keepStandardChromosomes(txDb)
    regionRanges <- rowRanges(experiments(object)[[experimentName]])
    beforeCols <- colnames(mcols(regionRanges))
    peakSeqLevs <- seqlevelsStyle(regionRanges)
    if( any(peakSeqLevs %in% "NCBI") )
        peakSeqLevs <- "Ensembl"
    seqlevelsStyle(bsGenome) <- peakSeqLevs
    seqlevelsStyle(txDb) <- peakSeqLevs
    ## subset txdb ##
    message("Creating txdb with active genes...")
    if( is.null(activeGeneList) ){
        message("No vector of active genes was specified, active genes are calculated based on peaks overlapping promoters.")
        promoters_gr <- promoters(genes(txDb, columns = columns(txDb)), use.names = TRUE)
        activeGeneList <- annotatePeakInBatch(
            regionRanges,
            AnnotationData = promoters_gr,
            output = "overlapping",
            select = "first" )$feature
        activeGeneList <- as.vector(na.omit(activeGeneList))
    }
    txDbActive <- subsetTxDbByGenes(txDb, activeGeneList)
    ## make annotation queries with diff params ##
    message("Annotating regions...")
    annotatePeakParams <- data.frame(
        annoName=c("all", "TSS", "AllFlank", "Active"),
        txDbName=c("txDb", "txDb", "txDb", "txDbActive"),
        overlap=c("all", "TSS", "all", "all"),
        frankDistance = c(5000, 5000, flankDistance, 5000),
        addFlankGeneInfo=c(FALSE, FALSE, TRUE, FALSE) )
    geneCentrAnnos <- lapply( seq_len(nrow(annotatePeakParams)), function(i){
        annotatePeak(
            regionRanges,
            TxDb=get(annotatePeakParams[i,"txDbName"]),
            tssRegion = c(-3000, 3000),
            verbose = FALSE,
            addFlankGeneInfo=annotatePeakParams[i,"addFlankGeneInfo"],
            annoDb = orgDb$packageName,
            overlap = annotatePeakParams[i,"overlap"] )
    } )
    names(geneCentrAnnos) <- annotatePeakParams$annoName
    ## Organize output of annotation queries ##
    resultRanges <- geneCentrAnnos[["all"]]@anno
    ## resultRanges <- regionRanges
    outParams <- data.frame(
        annoName=rep(c("TSS", "AllFlank", "Active"), c(2, 3, 3)),
        inName=c("SYMBOL", "geneId", "flank_txIds",
                 "flank_gene_distances",
                 "flank_geneIds", "geneId", "SYMBOL", "distanceToTSS"),
        outName=c("nearest_TSS", "nearest_TSSIds", "flank_txIds",
                  "flank_geneDistance", "flank_geneIds",
                  "active_geneIds", "active_geneSymb", "active_geneDistance") )
    for( i in seq_len(nrow(outParams)) ){
##        mcols(resultRanges)[[outParams$outName[i]]] <-
        ##            mcols(geneCentrAnnos[[outParams$annoName[i]]]@anno)[[outParams$inName[i]]]
        outRanges <- geneCentrAnnos[[outParams$annoName[i]]]@anno
        ovlR <- findOverlaps( resultRanges, outRanges, type="equal" )
        mcols(resultRanges)[[outParams$outName[i]]] <- NA
        if( outParams$inName[i] == "SYMBOL" & is.null(mcols(outRanges)$SYMBOL) ){
            mcols(outRanges)$SYMBOL <- mcols(outRanges)$geneId
        }
        mcols(resultRanges)[[outParams$outName[i]]][queryHits(ovlR)] <-
            mcols(outRanges)[[outParams$inName[i]]][subjectHits(ovlR)]
    }
    resultRanges$simple_annotation <- gsub(" *\\(.*?\\) *", "", resultRanges$annotation)
    if( is.null( chromHmmAnnotation ) ){
      if(download == TRUE){
        chromHmmAnnotation <- selectChromHMMData( object )
      }
    }
    if( !is.null( chromHmmAnnotation ) ){
        chromHMMData <- import( chromHmmAnnotation, format="bed" )
        names(chromHMMData) <- sprintf("annFeature%0.9d", seq_len(length(chromHMMData)))
        seqlevelsStyle(chromHMMData) <- peakSeqLevs
        chromCentrAnno <- annotatePeakInBatch(
            regionRanges,
            AnnotationData = chromHMMData,
            output = "overlapping", select = "first")
        resultRanges$chromHMM_annotation <- NA
        resultRanges$chromHMM_annotation[!is.na(chromCentrAnno$feature)] <-
            chromHMMData[na.omit(chromCentrAnno$feature)]$name
        resultRanges$chromHMM_annotation_simple <-
            simplifyChromHMMLabs( resultRanges$chromHMM_annotation )
    }
    ## make sure order is preserved ##
    ovl <- findOverlaps(regionRanges, resultRanges, type="equal")
    stopifnot(all(queryHits(ovl) == subjectHits(ovl)))
    ## annotate metadata ##
    afterCols <- colnames(mcols(resultRanges))
    if(is.null(elementMetadata(mcols(resultRanges))$type)){
        elementMetadata(mcols(resultRanges))$type <- NA
    }
    elementMetadata(mcols(resultRanges))$type[!afterCols %in% beforeCols] <- "Chraw2 annotation"
    rowRanges(experiments(object)[[experimentName]]) <- resultRanges
    object <- createUCSCLinks( object, experimentName )
    object
}
