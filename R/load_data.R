#' Function to rewrite paths in example datasets
#'
#' @description
#' This function rewrites the paths from the bam, bigwig, peak and qc files
#' available as test datasests for the package, allowing for the package to
#' find all the necessary data to run test code.
#'
#' @param ce A ChrawExperiment object
#' @export

 rewrite_paths <- function(ce) {
  inst_dir <- system.file("files/", package="chraw")
  ce$bamFile <- gsub("inst/files", inst_dir, ce$bamFile)
  ce$bigwigFile <- gsub("inst/files", inst_dir, ce$bigwigFile)
  ce$peakFile <- gsub("inst/files", inst_dir, ce$peakFile)
  ce$qcFile <- gsub("inst/files", inst_dir, ce$qcFile)
  return(ce)
 }
