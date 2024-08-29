#' @importMethodsFrom BiocGenerics mean
getPlottingData <- function( bigwigFile, plottingRegion, resolution=100, smoothIter=1 ){
  smoothWindowWidth <- resolution * 3
  plottingRegion <- resize(plottingRegion, width(plottingRegion)*1.5,fix="center")
  plottingRegion <- resize(plottingRegion, width(plottingRegion) + (smoothWindowWidth - (width(plottingRegion) %% smoothWindowWidth)), fix="center")
  plottingRegionTiles <- tile(plottingRegion, width=resolution)[[1]]
  smoothingTiles <- slidingWindows(plottingRegion, width=smoothWindowWidth, step=resolution)[[1]]
  covData <- import( bigwigFile, which=plottingRegion )
  seqlevels(covData) <- seqlevelsInUse(covData)
  binnedCov <- binnedAverage(plottingRegionTiles, coverage( covData, weight=covData$score ), varname="cov")
  while( smoothIter > 0 ){
    ovl <- findOverlaps(binnedCov, smoothingTiles)
    smoothingTiles$cov <- mean(NumericList(split(binnedCov$cov[queryHits(ovl)], subjectHits(ovl))))
    smoothedTiles <- resize( smoothingTiles, width=resolution, fix="center" )
    ovl2 <- findOverlaps( binnedCov, smoothedTiles, type="equal" )
    binnedCov$cov[queryHits(ovl2)] <- smoothedTiles$cov[subjectHits(ovl2)]
    binnedCov$cov[1] <- binnedCov$cov[2]
    binnedCov$cov[length(binnedCov)-1] <- binnedCov$cov[length(binnedCov)]
    smoothIter <- smoothIter - 1
  }
  binnedCov
}

#' @title Constructor of coverage data
#' @description This function creates a list of \code{GRanges} objects with the average coverage from
#' samples of a ChrawExperiment object. **NOTE**: samples names should end with replicate information
#' separated by a "_" to create the correct replicate groups.
#' @param object A ChrawExperiment object.
#' @param includeSamples Vector of strings with sample names to be included in coverage calculation.
#' @param plottingRegion A GRanges object with the region to be visualized.
#' @param resolution Numerical value for the resolution of the plot. Default value is 100.
#' @param smoothIter Number of iterations to smooth the plot. Default value is 2.
#'
#' @examples
#' data(ce_examples)
#' ce_examples <- rewrite_paths(ce_examples)
#'
#' includeSamples <- rownames(colData(ce_examples))[1:3]
#' plottingRegion <- GRanges("chr6:67090000-67170000")
#' covData <- importAndAverage(ce_examples,includeSamples,plottingRegion )
#'
#' @returns Dataframe of coverage data
#'
#' @export
importAndAverage <- function( object, includeSamples, plottingRegion, resolution=100, smoothIter=2 ){
  if(!is(object,"ChrawExperiment")){
    stop("object` must be a ChrawExperiment object")
  }
  possibleSamples <- rownames(colData(object))

  if(!is.null(includeSamples)){
    if( !all(includeSamples %in% rownames(colData(object)))) {
      stop(sprintf("The variables '%s' not found in the ChrawExperiment object. The `includeSamples=` parameter must be samples from `colData(object)`.", includeSamples))
    }
  }
  else{
    includeSamples <- possibleSamples
  }
  if(!is(plottingRegion,"GRanges")){
    stop("plottingRegion must be a GRanges object")
  }

  bigwigFiles <- colData(object)[rownames(colData(object)) %in% includeSamples,"bigwigFile"]
  covData <- lapply( bigwigFiles, function(xx){
    covData <- getPlottingData( xx, plottingRegion, resolution=resolution, smoothIter=smoothIter )
    covData <- as.data.frame(covData)
    covData
  })
  names(covData) <- colData(object)[rownames(colData(object)) %in% includeSamples,"sampleName"]
  for( i in names(covData) ){
    covData[[i]]$sample_alias <- i
  }
  covData <- do.call(rbind, covData)
  #Replicate group finds the last underline character and removes the information after it
  covData$replicate_group <- gsub("_[^_]*$", "", covData$sample_alias)

  covData <- aggregate(cov ~ seqnames + start + end + replicate_group, data = covData, FUN = mean)
  covData
}

#' @title Plot to visualize Genome Browser
#' @description This function creates the plot to visualize coverage data in a specific plotting region.
#' It takes as input the return value of the [importAndAverage] function.
#' @param covData Coverage data resulted from importAndAverage function.
#' @param replicate_group String containing the name of the replicate group to be plotted.
#' @param plottingRegion GRanges object of the plotting region to be visualized.
#' @param ylim Value of the limit in the y-axis. Default value is \code{NULL}.
#' @param cl String with the color of the bins. Default value is "black".
#' @param yexpand Logical. If TRUE, the default, adds a small expansion factor to the limits to ensure that data and axes don't overlap. If FALSE, limits are taken exactly from the data or ylim.
#'
#' @examples
#' data(ce_examples)
#' ce_examples <- rewrite_paths(ce_examples)
#'
#' includeSamples <- rownames(colData(ce_examples))[1:3]
#' plottingRegion <- GRanges("chr6:67090000-67170000")
#' covData <- importAndAverage(ce_examples,includeSamples, plottingRegion )
#' plot <- plotCovFromDF(covData,replicate_group="H3K27ac_NR1I3_CTRL_3h",plottingRegion, cl="red")
#'
#' @returns ggplot object
#'
#' @export
plotCovFromDF <- function( covData, replicate_group, plottingRegion, ylim=NULL, cl="black", yexpand=TRUE ){
  if(!is.null(covData)){
    if(!is(covData,"data.frame")){
      stop("covData must be a dataframe, resulted from importAndAverage function")
    }
  }
  else{
    stop("covData value is missing")
  }
  if(!is(plottingRegion,"GRanges")){
    stop("plottingRegion must be a GRanges object")
  }
  if(!is.null(replicate_group)){
    if(!(replicate_group %in% covData$replicate_group)){
      stop("replicate_group must be a string containing a valid replicate group from covData")
    }
  }
  else{
      stop("replicate_group value is missing")
  }
  covPlot <- covData[covData$replicate_group %in% replicate_group, ] %>%
    ggplot( ) +
    geom_rect(aes(xmin = start-1, xmax = end+1, ymax = cov, ymin=0 ), fill=cl) +
    scale_x_reverse( limits=c(end(plottingRegion), start(plottingRegion)), expand = c(0, 0)) +
    coord_cartesian( ylim=ylim, expand=yexpand ) +
    theme_cowplot() +
    labs(y=replicate_group)
  covPlot
}
