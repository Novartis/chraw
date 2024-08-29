getDedupAlignments <- function(x){
    stopifnot(all(file.exists( x )))
    x <- gsub(".bam$", ".nodup.bam", x)
    stopifnot(all(file.exists(x)))
    x
}

getIdrPeakFile <- function( x ){
    fl <- file.path( dirname(dirname( x )), "idr_reproducibility", "idr.optimal_peak.narrowPeak.gz" )
    ifelse( file.exists(fl), fl, NA )
}

getOverlapPeakFile <- function( x ){
    fl <- file.path( dirname(dirname( x )), "overlap_reproducibility", "overlap.optimal_peak.narrowPeak.gz" )
    ifelse( file.exists(fl), fl, NA )
}

#' @importFrom jsonlite fromJSON
readQcJson <- function( basepaths, pipeline, subDirName, qcFiles ){
  if( pipeline %in% c("chip-seq", "atac-seq") ){
    basepaths <- na.omit( basepaths )
    qcFiles <- unique(file.path( basepaths, "qc", "qc.json" ))
    names(qcFiles) <- basename(dirname(dirname(qcFiles)))
  }else if( pipeline != 'public' ){
    basepaths <- na.omit( basepaths )
    qcFiles <- unique(file.path( basepaths, "qc", subDirName, "qc.json" ))
    names(qcFiles) <- basename(dirname(qcFiles))
  }
  existsFlag <- file.exists(qcFiles)
  if( !all(existsFlag) ){
    qcFiles <- file.path( dirname( qcFiles ), sprintf("%s.qc.json", basename(dirname(qcFiles))))
    names(qcFiles) <- basename(dirname(qcFiles))
    existsFlag <- file.exists(qcFiles)
  }
  if( !all(existsFlag) ){
    missingFiles <- qcFiles[!existsFlag]
    plural <- ifelse(length(missingFiles) > 1, "s", "")
    tns <- ifelse( plural == "s", "are", "is" )
    stop(sprintf("The following quality control json file%s %s missing from the cromwell output%s:\n\t%s",
                 plural, tns, plural, paste(missingFiles, collapse="\n\t")))
  }
  qcData <- lapply(qcFiles, function(x){
    jsonlite::fromJSON(x)
  } )
  names(qcData) <- names(qcFiles)
  qcData
}

readQcJsonFromCE <- function( ce ){
  if( ce@pipeline == 'public' ){
    qcFiles <- unique(colData(ce)$qcFile)
    names(qcFiles) <- unique(colData(ce)$condition)
    readQcJson(qcFiles = qcFiles,
               pipeline = ce@pipeline)
  }else{
    readQcJson(basepaths = colData(ce)$cromwellBaseOutput,
               pipeline = ce@pipeline,
               subDirName = colData(ce)$`Replicate Group` )
  }
}

getCromwellGenome <- function( basepaths, pipeline, subDirName ){
    qcData <- readQcJson( basepaths, pipeline, subDirName )
    refGenomeVector <- vapply( qcData, function(x){x[["general"]]$genome}, character(1) )
    refGenomeVector
}

#' Import peaks from *.narrowPeak file into GRange object
#'
#' @param narrowPeak.file File path to the narrow peaks bed file.
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#'
#' @return GRange object
#' @export
importNarrowPeaksFromFile <- function(narrowPeak.file){
    peaks <- read.table(narrowPeak.file, sep = "\t")
    colnames(peaks) <- c("seqname", "start", "end", "name",
                         "score", "strand", "signalValue",
                         "pValue", "qValue", "peak")
    # Change column types
    peaks$seqname <- as.character(peaks$seqname)
    peaks$start <- as.numeric(peaks$start)
    peaks$end <- as.numeric(peaks$end)
    peaks$name <- as.character(peaks$name)
    peaks$score <- as.numeric(peaks$score)
    peaks$strand <- as.character(peaks$strand)
    peaks$signalValue <- as.numeric(peaks$signalValue)
    peaks$pValue <- as.numeric(peaks$pValue)
    peaks$qValue <- as.numeric(peaks$qValue)
    peaks$peak <- as.numeric(peaks$peak)

    makeGRangesFromDataFrame( peaks, keep.extra.columns = TRUE )
}

#' Import peaks from a ChrawExperiment object
#'
#' @description This function inputs a ChrawExperiment and imports
#' the peak for the specified samples. The output is a GRangesList
#' by default.
#'
#' @param object A ChrawExperiment object.
#' @param includeSamples A character vector specifying the samples to be included in the new experiment. These names should match with `rownames(colData(object))`.
#' @param peakType A character string. Valid values are 'perSample' to import the peaks called per sample, "idr" to import the IDR peak calls and "overlap" to
#' @param merge A logical flag indicating whether peaks should be merged.
#' @return A ChrawExperiment object with a new experiment added.
#' @importFrom Rsubread featureCounts
#' @import BiocParallel
#' @importFrom random randomStrings
#' @export
importNarrowPeaks <- function( object, includeSamples=rownames(colData(object)), peakType="perSample", merge=FALSE ){
    includeFlag <- includeSamples %in% rownames( colData( object ) )
    if( !all( includeFlag ) ){
        stop(sprintf("The following sample(s) could not be found:\n\t%s\n       Please make sure that the samples match with `rownames(colData(object))`.",
                     paste(includeSamples[!includeFlag], collapse="\n\t" )))
    }
    ## check and locate peak files depending on the peakType ##
    if( !peakType %in% c("perSample", "idr", "overlap") ){
        stop("Invalid 'peakType' parameter. 'peakType' must be either 'perSample', 'idr' or 'overlap'")
    }
    peakFiles <- colData(object)$peakFile
    names(peakFiles) <- rownames(colData(object))
    locatorFun <- switch(
        peakType,
        "perSample" = function(x){x},
        "idr" = getIdrPeakFile,
        "overlap" = getOverlapPeakFile )
    peakFiles <- locatorFun( peakFiles )
    names(peakFiles) <- rownames(colData(object))
    peakFiles <- peakFiles[includeSamples]
    existsFlag <- file.exists( peakFiles )
    if( !all( existsFlag ) ){
        warning(sprintf("The peak files for the following samples could not be found in the file system and will be ignored:\n\t%s\n",
                        paste(names(peakFiles)[!existsFlag], collapse="\n\t")))
        peakFiles <- peakFiles[existsFlag]
    }
    if( length(peakFiles) == 0 ){
        stop("No peak files were found")
    }
    peaksPerFile <- lapply( peakFiles, importNarrowPeaksFromFile )
    names(peaksPerFile) <- names(peakFiles)
    if( merge ){
        peaksPerFile <- reduce(Reduce(c, peaksPerFile))
    }else{
        peaksPerFile <- GRangesList(peaksPerFile)
        names(peaksPerFile) <- names(peakFiles)
    }
    peaksPerFile
}
