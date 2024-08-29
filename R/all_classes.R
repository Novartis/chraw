#' @title ChrawExperiment class
#' @description This class is an extension of the MultiAssayExperiment class,
#' containing also information on the reference genome from the experiment
#' and the pipeline used to create the experiments stored in the objects.
#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @export
setClass(
    "ChrawExperiment",
    contains = "MultiAssayExperiment",
    slots=list( referenceGenome="character", pipeline="character" ) )

setValidity(
    "ChrawExperiment",
    function(object) {
        supportedGenomes <- c("hg38", "hg19", "mm10", "mm9", "rn6", "canFam3")
        requiredColumns <- c("bamFile", "bigwigFile", "peakFile")
        if(!object@referenceGenome %in% supportedGenomes){
            return(sprintf("The reference genome, '%s', is currently not supported", object@referenceGenome))
        }
        supportedPipelines <- c("public", "chromatinalign", "chip-seq", "atac-seq")
        if(!object@pipeline %in% supportedPipelines){
            return(sprintf("The pipeline, '%s', is currently not supported", object@pipeline))
        }
        if( length(object@referenceGenome) != 1 ){
            return("The reference genome slot must be a character vector of length 1")
        }
        rFlag <- requiredColumns %in% colnames(colData(object))
        if( !all(rFlag) )
            return(sprintf("The following column(s) are missing in the colData:\n\t%s", paste(requiredColumns[!rFlag], collapse="\n\t")))
        if( !"ExistingSampleFlag" %in% names(experiments(object) ) )
            return("An experiment 'ExistingSampleFlag' indicating the presence of a sample is missing")
        TRUE
    }
)

#' @title Constructor of ChrawExperiments from a dataframe
#' @description This function builds a ChrawExperiment object from a pre-defined
#' dataframe containing following information:
#' 1) sampleName - Sample names.
#' 2) bamFile - Path to alignment file in .bam format.
#' 3) bigwigFile - Path to signal tracks in .bigwig format.
#' 4) peakFile - Path to peaks in either bed, narrowPeaks or broadPeak format.
#' 5) qcFile - Path to qc file in json format (output from ENCODE pipeline).
#' 6-) Any additional sample information.
#' The ChrawExperiment object is an extension of the MultiAssayExperiment
#' object specifically designed for the analysis of chromatin experiments.
#' The ChrawExperiment object aims at being the central data structure for
#' downstream analyses.
#'
#' @param sampleDF Dataframe with required columns
#' @param genome Reference genome can be: mm9, mm10, hg19, hg38, rn6.
#' @param minimalMetadata Logical, indicating whether all metadata from the sample sheet file should be kept or only a set of essential columns.

#' @examples
#' sampleFile <- system.file("files", "sampleSheet.txt", package="chraw", mustWork=TRUE)
#' sampleDF <- as.data.frame(read.csv(sampleFile, sep = "\t", header=TRUE))
#'
#' inst_dir <- system.file("files/", package="chraw")
#' sampleDF$bamFile <- gsub("inst/files", inst_dir, sampleDF$bamFile)
#' sampleDF$bigwigFile <- gsub("inst/files", inst_dir, sampleDF$bigwigFile)
#' sampleDF$peakFile <- gsub("inst/files", inst_dir, sampleDF$peakFile)
#' sampleDF$qcFile <- gsub("inst/files", inst_dir, sampleDF$qcFile)
#'
#' ce <- ChrawExperimentFromDataFrame( sampleDF , genome = 'rn6')
#'
#' @return A ChrawExperiment object.
#' @export
ChrawExperimentFromDataFrame <- function(sampleDF, genome, minimalMetadata=TRUE){
  sampleData <- readAndValidateDataFrame(sampleDF, minimalMetadata=minimalMetadata )
  pipeline <- 'public'
  # Check if files exist
  fileCols <- colnames(sampleData)[grep(colnames(sampleData), pattern = 'File')]
  for( fileCol in fileCols ) {
    fileExists <- !sapply(sampleData[,fileCol], file.exists)
    if( sum(fileExists) > 0 ) {
      stop(sprintf("Could not find following %s files: %s",
                   fileCol, sampleData[fileExists,fileCol]))
    }
  }
  # Create ChrawExperiment
  mat <- matrix( rep(TRUE, times=nrow(sampleData)), nrow=1L )
  colnames(mat) <- sampleData$sampleName
  rownames(sampleData) <- colnames(mat)
  cr <- new("ChrawExperiment",
            MultiAssayExperiment(
              list(ExistingSampleFlag=mat),
              colData=DataFrame(sampleData, check.names=FALSE) ),
            referenceGenome=genome,
            pipeline=pipeline )
  rownames(colData(cr)) <- colData(cr)$sampleName
  validObject(cr)
  cr
}


#' @title Constructor of ChrawExperiments from a sample sheet
#' @description This function builds a ChrawExperiment object from a pre-define
#' sample sheet containing following information:
#' 1) sampleName - Sample names.
#' 2) bamFile - Path to alignment file in .bam format.
#' 3) bigwigFile - Path to signal tracks in .bigwig format.
#' 4) peakFile - Path to peaks in either bed, narrowPeaks or broadPeak format.
#' 5-) Any additional sample information.
#' The ChrawExperiment object is an extension of the MultiAssayExperiment
#' object specifically designed for the analysis of chromatin experiments.
#' The ChrawExperiment object aims at being the central data structure for
#' downstream analyses.
#'
#' @param sampleFile Path to the sample sheet file.
#' @param genome Reference genome can be: mm9, mm10, hg19, hg38, rn6.
#' @param minimalMetadata Logical, indicating whether all metadata from the sample sheet file should be kept or only a set of essential columns.

#' @examples
#' #sampleFile <- system.file("files", "sampleSheet.txt", package="chraw", mustWork=TRUE)
#' #ce <- ChrawExperimentFromSampleFile( sampleFile , genome = 'rn6')
#'
#' @return A ChrawExperiment object.
#' @export
ChrawExperimentFromSampleFile <- function( sampleFile, genome, minimalMetadata=TRUE ){
    sampleData <- readAndValidateSampleSheet( sampleFile, minimalMetadata=minimalMetadata )
    pipeline <- 'public'
    # Check if files exist
    fileCols <- colnames(sampleData)[grep(colnames(sampleData), pattern = 'File')]
    for( fileCol in fileCols ) {
        fileExists <- !sapply(sampleData[,fileCol], file.exists)
        if( sum(fileExists) > 0 ) {
            stop(sprintf("Could not find following %s files: %s",
                         fileCol, sampleData[fileExists,fileCol]))
        }
    }

    # Create ChrawExperiment
    mat <- matrix( rep(TRUE, times=nrow(sampleData)), nrow=1L )
    colnames(mat) <- sampleData$sampleName
    rownames(sampleData) <- colnames(mat)
    cr <- new("ChrawExperiment",
              MultiAssayExperiment(
                  list(ExistingSampleFlag=mat),
                  colData=DataFrame(sampleData, check.names=FALSE) ),
              referenceGenome=genome,
              pipeline=pipeline )
    rownames(colData(cr)) <- colData(cr)$sampleName
    validObject(cr)
    cr
}
