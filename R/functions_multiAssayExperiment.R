#' Add fragment length distribution per sample to ChrawExperiment.
#'
#' @description This function calculates the length distribution
#' of the sequenced fragments and adds these information as an experiment
#' to the ChrawExperiment object.
#'
#' @param object A ChrawExperiment object.
#' @param maxFragLength A numeric value specifying the maximum fragment length to consider.
#' @param BPPARAM A BiocParallel instance for parallelization.
#' @param ... Further options passed to `csaw::getPESizes()`
#'
#' @import MultiAssayExperiment
#' @importFrom BiocParallel bplapply
#'
#' @return ChrawExperiment object with "FragLengthDist" experiment
#' @export
#'
addFragmentLengthDist <- function( object, maxFragLength = 750,
                                   BPPARAM = SerialParam(), ... ) {
  if( !is( object, "ChrawExperiment" ) )
    stop("Parameter 'object' needs to be a ChrawExperiment object")

  validObject(object)

  # Check if *.bams are paired-end
  qcData <- readQcJsonFromCE( object )
  isPaired <- sapply(qcData, function(x) x[['general']][['seq_endedness']][[1]])
  if( sum(unlist(isPaired)) != length(isPaired) )
    stop("Can't compute fragment length distribution from single end bam files")

  ## Get input *.bams
  bams <- colData(object)$bamFile
  names(bams) <- rownames(colData(object))

  ## Compute fragment length for all experiments
  fragLen.dists <- bplapply(
    bams,
    function(x, ...){
      computeFragLengthDist(x, ...)
    }, BPPARAM = BPPARAM , ... )

  # fragLen.dists <- bplapply(
  #   bams,
  #   function(x){
  #     computeFragLengthDist(x, param=csaw::readParam(pe = "both", restrict="chr1"))
  #   }, BPPARAM = MulticoreParam(2) )
  # names(fragLen.dists) <- names(bams)

  ## Compute fragment length frequency matrix from tidy fragment length dists
  fragLen.dists <- do.call(rbind, fragLen.dists)
  fragLen.dists$bamReads <- names(bams)[match( fragLen.dists$bamReads, bams )]
  fragLen.dists <- fragLen.dists[(fragLen.dists$fraglen <= maxFragLength),]
  fragLen.dists <- split( fragLen.dists, fragLen.dists$bamReads )
  fragLen.dists <- fragLen.dists[names(bams)]
  mat <- matrix( 0, ncol=length(bams), nrow=maxFragLength )
  colnames(mat) <- names(bams)
  for( i in names(fragLen.dists) ){
    mat[fragLen.dists[[i]]$fraglen,i] <- fragLen.dists[[i]]$freq
  }

  ## ass assay to experiment
  addExperiment( object, mat, "FragLengthDist" )
}

#' Adds a count experiment to a ChrawExperiment object
#'
#' @description This function inputs a ChrawExperiment and counts the number
#' of read fragments overlapping with a set of specified genomic ranges. The read
#' counting can be done either for all samples or only a subset of them. The results
#' are formatted as a SummarizedExperiment object and saved as part of the
#' ChrawExperiment object.
#'
#' @param object A ChrawExperiment object.
#' @param regions A named GenomicRegions with the genomic coordinates of the regions to count.
#' @param name A string specifying the name of the new experiment.
#' @param includeSamples A character vector specifying the samples to be included in the new experiment. These names should match with `rownames(colData(object))`.
#' @param dedup Logical indicating whether to use the bam files without PCR duplicates or not.
#' @param isPairedEnd Logical indicating whether reads are paired-end or not.
#' @param normalizeCPM Logical indicating whether counts should be normalized to CPMs or not.
#' @param BPPARAM A BiocParallel instance.
#' @param ... additional arguments passed to `featureCounts()`
#' @return A ChrawExperiment object with a new experiment added.
#' @importFrom Rsubread featureCounts
#' @import BiocParallel
#' @importFrom random randomStrings
#' @export
addCountExperiment <- function( object, regions, name, includeSamples=rownames(colData(object)),
                                dedup=FALSE, isPairedEnd=TRUE, normalizeCPM=TRUE, BPPARAM=SerialParam(), ... ){
  validObject(object)
  bamFiles <- colData(object)$bamFile
  names(bamFiles) <- rownames(colData(object))
  includeFlag <- includeSamples %in% names(bamFiles)
  if( !all( includeFlag ) ){
    stop(sprintf("The following sample(s) could not be found:\n\t%s\n       Please make sure that the samples match with `rownames(colData(object))`.",
                 paste(includeSamples[!includeFlag], collapse="\n\t" )))
  }
  bamFiles <- bamFiles[includeSamples]
  if( any( is.na(bamFiles) ) ){
    if(all(is.na(bamFiles))){
      stop("All bam file paths are set to NAs.")
    }
    warning(sprintf("Skipping samples having NA values as bam file paths: %s",
                    paste(names(bamFiles)[is.na(bamFiles)], collapse=",")))
    bamFiles <- bamFiles[!is.na(bamFiles)]
  }
  if( dedup ){
    bamFiles <- getDedupAlignments(bamFiles)
  }
  if( !all(file.exists(bamFiles) ))
    stop("At least one of the specified bam files does not exist")
  if(any(duplicated( basename(bamFiles) )))
    stop("Bam file names are not unique")
  if(is.null(names(regions)))
    stop("The 'regions' parameter needs to be a named GRanges object")
  if(any(duplicated(names(regions))))
    stop("Names of 'regions' must be unique")
  if(any(duplicated(regions)))
    stop("Each element of 'regions' must be unique")
  if(!is( name, "character" ) )
    stop("Parameter 'name' must be a character vector")
  ### count reads using featureCounts ###
  binDF <- data.frame(
    GeneID = names(regions),
    Chr = as.character( seqnames( regions ) ),
    Start = start( regions ),
    End = end( regions ),
    Strand = "+" )
  countData <- bplapply( bamFiles, function( x, binDF, ... ){
    cnts <- Rsubread::featureCounts(
      x,
      annot.ext = binDF,
      isPairedEnd = isPairedEnd,
      countMultiMappingReads = FALSE,
      allowMultiOverlap = TRUE,
      nthreads = 1, ... )
    cnts
  }, BPPARAM=BPPARAM, binDF=binDF, ... )
  ## reformat output into a SummarizedExperiment object
  countStats <- do.call(cbind, lapply( countData, function(x){x$stat[,2L,drop=FALSE]} ))
  colnames(countStats) <- names(countData)
  countStats <- t(countStats)
  colnames(countStats) <- countData[[1L]]$stat[,"Status"]
  countMatrix <- do.call(cbind, lapply( countData, "[[", "counts" ))
  colnames( countMatrix ) <- names( countData )
  expList <- list('counts' = countMatrix)
  ## Normalize counts to CPM
  if(normalizeCPM) {
    libSize <- rowSums(as.matrix(countStats))
    cpmMatrix <- sapply(names(libSize), function(sample) {
      10^6*(countMatrix[,sample]/libSize[sample])
    })
    expList[['cpms']] <- cpmMatrix
  }
  ## Combine in SummarizedExperiment
  se <- SummarizedExperiment( assays=expList, colData=DataFrame( countStats ),
                              rowRanges=regions )
  colnames(se) <- names(countData)
  ## combine to ChrawExperiment ##
  object <- addExperiment( object, se, name )
  object
}

addExperiment <- function( object, se, name, sampleMap=NULL ){
  se <- SimpleList( se )
  if( name %in% names(experiments(object)) ){
    message(sprintf("An experiment named %s was found. Replacing this with the new experiment...\n", name))
    tmpname <- as.vector(randomStrings(1, 15))
    names(se) <- tmpname
    object <- c(object, se, sampleMap=sampleMap)
    experiments(object)[[name]] <- NULL
    names(object)[names(object) %in% tmpname] <- name
  }else{
    names(se) <- name
    object <- c(object, se, sampleMap=sampleMap)
  }
  object
}

#' Add an RNA-seq experiment to a ChrawExperiment object
#'
#' @description This function inputs a ChrawExperiment object
#' and either a RangedSummarizedExperiment or a DESeqDataSet object
#' containing RNA-seq count data. The function integrates the
#' RNA-seq count data into the ChrawExperiment container. It also
#' maps ENTREZ, ENSEMBL and SYMBOL identifiers to enable
#' integration with the chromatin datasets.
#'
#' @param object A ChrawExperiment object.
#' @param se Either a 'RangedSummarizedExperiment' or a 'DESeqDataSet' object.
#' @param name A string specifying the name of the new experiment.
#' @param identifierType A character string specifying the identifier type. Can be 'ENSEMBL', 'ENTREZID' or 'SYMBOL'.
#' @return A ChrawExperiment object with an RNA-seq experiment added.
#'
#' @importFrom ensembldb select
#'
#' @export
addRNASeqExperiment <- function( object, se, name=NULL, identifierType="ENSEMBL" ){
  if( !is( object, "ChrawExperiment") )
    stop("Parameter 'object' must be a ChrawExperiment object")
  if( !(is( se, "DESeqDataSet") | is( se, "RangedSummarizedExperiment" )) ){
    stop("Parameter 'se' must be either a DESeqDataSet or a RangedSummarizedExperiment object")
  }
  if( !is( name, "character") )
    stop("The parameter 'name' must be a character vector")
  if( is( se, "DESeqDataSet" ) ){
    rowNms <- rownames(se)
    se <- as(se, "RangedSummarizedExperiment")
    rownames( se ) <- rowNms
  }
  if( is.null( rownames(se) ) )
    stop("The rownames of the `se` object are not specified. The `rownames` of the `se` object must be named by gene identifiers.")
  nameTypes <-c("ENTREZID", "SYMBOL", "ENSEMBL")
  if( !identifierType %in% nameTypes )
    stop("Parameter 'identifierType' needs to be either 'ENTREZID', 'SYMBOL' or 'ENSEMBL'")
  if( !is(rowRanges(se), "GRanges") )
    stop("The object 'se' must have the gene coordinates as GRanges in the rowRanges field.")
  orgdb <- selectOrgDb(object)
  geneMaps <- select(orgdb, keys=rownames(se), columns=nameTypes, keytype=identifierType)
  geneMapsDF <- DataFrame(geneMaps[[identifierType]])
  colnames(geneMapsDF) <- identifierType
  rownames(geneMapsDF) <- geneMapsDF[[identifierType]]
  for( x in nameTypes[!nameTypes %in% identifierType] ){
    geneMapsSub <- na.omit(geneMaps[,c(x, identifierType)])
    sp <- split( geneMapsSub[[x]], geneMapsSub[[identifierType]] )
    geneMapsDF[[x]] <- NA
    geneMapsDF[names(sp),x] <- CharacterList(sp)
  }
  rowData(se) <- cbind( rowData(se), geneMapsDF[rownames(se),])
  se <- SimpleList(se)
  names(se) <- name
  se <- MultiAssayExperiment(
        experiments=ExperimentList(se),
        colData=colData(se[[name]]))
    refGen <- object@referenceGenome
    pipeline <- object@pipeline
    object <- c(object, se)
    attributes(experiments(object)[[name]])$metadata$isrnaseq <- TRUE
    new(
        "ChrawExperiment",
        object,
        referenceGenome=refGen,
        pipeline=pipeline )
}
