#' @importFrom S4Vectors SimpleList DataFrame
readAndValidateSampleSheet <- function( x, minimalMetadata=TRUE ){
    if( !file.exists( x ) )
        stop(paste("The sample sheet was not found. Are you specifying",
                   "correctly the path to the sample sheet file?"))
    sampleData <- read.delim( x, check.names=FALSE )
    requiredColumns <- c( "sampleName", "bamFile", "bigwigFile",
                         "peakFile","qcFile","condition")
    colFlag <- requiredColumns %in% colnames(sampleData)
    if( !all( colFlag ) ){
        stop(sprintf(paste("At least one of the following required columns",
                           "are missing from the sample sheet file: %s\n\t"),
                     paste(requiredColumns[!colFlag], sep="\n\t")))
    }
    ## nest rows with the same sample names, i.e. technical replicates
    if( any( duplicated( sampleData$sampleName ) ) ){
        stop(paste("Sample sheet contains multiple samples with the same sample",
                   "name. In case of technical replicates, please merge manually"))
    }
    if( minimalMetadata ){
        sampleData <- sampleData[,requiredColumns]
    }
    sampleData
}

#' @importFrom S4Vectors SimpleList DataFrame
readAndValidateDataFrame <- function(x, minimalMetadata=TRUE){
  sampleData <- data.frame(x)
  requiredColumns <- c( "sampleName", "bamFile", "bigwigFile",
                        "peakFile","qcFile","condition")
  colFlag <- requiredColumns %in% colnames(sampleData)
  if( !all( colFlag ) ){
    stop(sprintf(paste("At least one of the following required columns",
                       "are missing from the sample sheet file: %s\n\t"),
                 paste(requiredColumns[!colFlag], sep="\n\t")))
  }
  ## nest rows with the same sample names, i.e. technical replicates
  if( any( duplicated( sampleData$sampleName ) ) ){
    stop(paste("Sample sheet contains multiple samples with the same sample",
               "name. In case of technical replicates, please merge manually"))
  }
  if( minimalMetadata ){
    sampleData <- sampleData[,requiredColumns]
  }
  data.frame(sampleData)
}
