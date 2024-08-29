#' @title Peak to peak assignment function.
#' @description Assign peaks to peaks either within the same experiment or
#' across experiments using an assignment function such as `findOverlaps`,
#' `nearest`, etc.
#' @param object A ChrawExperiment object.
#' @param experimentNames Vector with two count experiment names, normally added using the `addCountExperiment` function.
#' @param filterAnnoCols Vector with one/two column names in `rowData` for peak filtering, e.g. simple_annotation added by the `annotateExperimentRegions` function.
#' @param filterValues  Vector with one/two values to bes used to filter the peaks, e.g. in case of simple_annotation it can be 'Promoter', 'Intron' etc.
#' @param assignFunc GenomicRanges function to compare peak sets. Add the moment supported functions include `findOverlaps`, `precede`, `follow` and `nearest`.
#'
#' @return A data.frame combining the rowData from experiment 1 with experiment 2 based on the peak assignment.
#'
#' @export
#'
#' @import GenomicRanges
#'
#' @examples
#' data(ce_examples)
#' ce_examples <- rewrite_paths(ce_examples)
#'
#' # Combine promoters with the nearest peak
#' assignPeaks2Peaks(object = ce_examples,
#'                   experimentNames = c('Peaks', 'Peaks'),
#'                   filterAnnoCols = c('simple_annotation', 'simple_annotation'),
#'                   filterValues = c('Promoter', 'Distal Intergenic'),
#'                   assignFunc = GenomicRanges::nearest)
#'
assignPeaks2Peaks <- function(object, experimentNames,
                              filterAnnoCols = NULL,
                              filterValues = NULL,
                              assignFunc = findOverlaps) {

  # Check if experimentNames are part of the obj
  if(sum(experimentNames %in% names(object)) != 2) {
    stop(sprintf('No experiment called %s in ChrawExperiment object. ',
                 experimentNames[!experimentNames %in% names(object)]))
  }

  # Check if experiments are RangedSummarizedExperiment
  exp.class <- sapply(experimentNames, function(exp) class(object[[exp]])[1])
  if(sum(exp.class == 'RangedSummarizedExperiment') != 2) {
    stop(sprintf('%s in ChrawExperiment object is not a RangedSummarizedExperiment. ',
                 experimentNames[exp.class != 'RangedSummarizedExperiment']))
  }

  # Get peaks for each experiment
  peaks <- lapply(experimentNames, function(expn) object[[expn]])

  # Filter peaks for annotation if defined
  if(!is.null(filterAnnoCols)) {
    peaks <- lapply(seq_len(length(filterAnnoCols)), function(ii) {
      if(!(filterAnnoCols[ii] %in% colnames(mcols(peaks[[ii]]))))
        stop(sprintf('%s is not existing in the rowData of the experiment %s.',
                     filterAnnoCols[ii], experimentNames[ii]))
      if(is.null(filterValues[ii]))
        stop(sprintf('Please specify a filterValue to subset the peaks for %s.',
                     filterAnnoCols[ii]))
      if(!(filterValues[ii] %in% mcols(peaks[[ii]])[,filterAnnoCols[ii]]))
        stop(sprintf('%s is not part of the filter colunmn %s.',
                    filterValues[ii], filterAnnoCols[ii]))
      peak <- peaks[[ii]]
      peak <- peak[mcols(peak)[,filterAnnoCols[ii]] == filterValues[ii],]
      return(peak)
    })
  }

  # Check if assignFunc is valid
  valid.funcs <- c('nearest', 'findOverlaps', 'precede', 'follow')
  if(!isS4(assignFunc)) {
    stop(paste0('The peak to peak assignment function is not a valid function. ',
                'Please use: ', paste(valid.funcs, collapse = ', ')))
  } else if(!assignFunc@generic[1] %in% valid.funcs) {
    stop(paste0('The peak to peak assignment function is not a valid function. ',
                'Please use: ', paste(valid.funcs, collapse = ', ')))
  }

  # Assign peaks to each other using different functions
  ov <- assignFunc(peaks[[1]], peaks[[2]])

  if(isS4(ov)) {
    cbind(mcols(peaks[[1]])[queryHits(ov),],
          mcols(peaks[[2]])[subjectHits(ov),])
  } else {
    cbind(mcols(peaks[[1]]), mcols(peaks[[2]])[ov,])
  }
}

#' @title Peak to Gene assignment function.
#' @description Assign peaks to genes either within the same experiment or
#' across experiments using an assignment function such as `findOverlaps`,
#' `nearest`, etc. Function should be used with the second experiment added
#' by `addRNASeqExperiment()` function.
#' @param object A ChrawExperiment object.
#' @param experimentNames Vector with two count experiment names, normally added using the `addCountExperiment` function.
#' @param filterAnnoCols Vector with one/two column names in `rowData` for peak filtering, e.g. simple_annotation added by the `annotateExperimentRegions` function.
#' @param filterValues  Vector with one/two values to bes used to filter the peaks, e.g. in case of simple_annotation it can be 'Promoter', 'Intron' etc.
#' @param contrasts Contrasts created from `testForDiffSignal` function for differential expression analysis.
#' @param geneAnnoCols Annotation column for genes.
#'
#' @return A data.frame combining the rowData from experiment 1 with experiment 2 based on the peak assignment.
#'
#' @export
#'
assignPeaks2GeneDE <-  function( object, experimentNames, contrasts, geneAnnoCols="nearest_TSS", filterAnnoCols=NULL, filterValues=NULL ){
    if( length(geneAnnoCols) != 1 | !is(geneAnnoCols, "character") )
        stop("If one of the experiments is an RNA-seq experiment, the 'geneAnnoCols' para must be a character vector of length 1.")
    if(sum(experimentNames %in% names(object)) != 2) {
        stop(sprintf('No experiment called %s in ChrawExperiment object. ',
                     experimentNames[!experimentNames %in% names(object)]))
    }
    exp1 <- experiments(object)[[experimentNames[1]]]
    exp2 <- experiments(object)[[experimentNames[2]]]
    if( is.null(metadata(exp2)$isrnaseq) )
        stop("The second experiment name must be an experiment added using the `addRNASeqExperiment()` function.",
             call.=FALSE)
    peakDat <- mcols(exp1)
    rnaseqDat <- mcols(exp2)
    if( (length(filterValues) != 1 | !is(filterValues, "character")) & !is.null(filterValues) )
        stop("If one of the experiments is an RNA-seq experiment, the 'filterValues' para must be a character vector of length 1.")
    if( !is.null(filterAnnoCols) ){
        if( (length(filterAnnoCols) != 1 | !is(filterAnnoCols, "character")) ){
            stop("If one of the experiments is an RNA-seq experiment, the 'filterAnnoCols' para must be a character vector of length 1.")
        }
        if(!filterAnnoCols %in% colnames(peakDat)){
            stop(sprintf("The column '%s' was not found in the rowData of the peak experiment.", filterAnnoCols))
        }
        peakDat <- peakDat[peakDat[[filterAnnoCols]] == filterValues,]
    }
    if(!geneAnnoCols %in% colnames(peakDat))
        stop(sprintf("The column '%s' was not found in the rowData of the peak experiment.", geneAnnoCols))
    ## Check if experimentNames are part of the obj
    if( nrow(peakDat) == 0 )
        stop("Nothing to plot. All peak data were filtered out by 'filterAnnoCols == filterValues'", call.=FALSE)
    numbMatch <- vapply(c("ENTREZID", "ENSEMBL", "SYMBOL"), function(x){
        sum(peakDat[[geneAnnoCols]] %in% rnaseqDat[[x]], na.rm=TRUE)
    }, numeric(1))
    numbMatch <- names(numbMatch)[which.max(numbMatch)]
    message(sprintf("Inferring identifier types for the column 'active_geneIds': %s", numbMatch))
    rnaMatch <- match(peakDat[[geneAnnoCols]], rnaseqDat[[numbMatch]])
    matchedDat <- cbind(peakDat, rnaseqDat[rnaMatch,])
    matchedDat
}


#' @title Scatterplot comparison of log2(FC) from differential analysis.

#' @description Assign peaks to peaks either within the same experiment or
#' across experiments using an assignment function such as `findOverlaps`,
#' `nearest`, etc.
#' @param object A ChrawExperiment object.
#' @param experimentNames Vector with two count experiment names, normally added using the `addCountExperiment` function.
#' @param contrasts Vector with two contrast names, used to compute diff peak statistics using the `testForDiffSignal` function.
#' @param filterAnnoCols Vector with one/two column names in `rowData` for peak filtering, e.g. simple_annotation added by the `annotateExperimentRegions` function.
#' @param filterValues  Vector with one/two values to bes used to filter the peaks, e.g. in case of simple_annotation it can be 'Promoter', 'Intron' etc.
#' @param assignFunc GenomicRanges function to compare peak sets. Add the moment supported functions include `findOverlaps`, `precede`, `follow` and `nearest`.
#' @param geneAnnoCols If one of the experiments is an RNA-seq dataset, the param `geneAnnoCols` must be a character vector that matches a column name of `mcols()` of the other experiment.
#' @param colorByFDR A string indicating how to color the points. Values can be either 'mean' (default) to color the point according to the mean q-value (FDR) across experiments, `exp1` and `exp2` to color by the q-value of individual experiments or 'none'.
#'
#' @return Scatterplot comparing the log2(FC) of two diff. contrasts/region sets and/or experiments.
#' @export
#'
#' @import ggplot2
#' @import cowplot
#' @import rlang
#'
#' @examples
#' data(ce_examples)
#' ce_examples <- rewrite_paths(ce_examples)
#'
#' # Perform diff testing for the scatterplot
#' contrasts <- list(contrast01 = c("condition", "agonist_3h", "CTRL_3h"),
#'                  contrast02 = c("condition", "agonist_27h", "CTRL_27h"))
#'
#' ce_examples <- testForDiffSignal(object = ce_examples,
#'                            experimentName = "Peaks",
#'                            design = ~condition,
#'                            contrasts = contrasts)
#'
#' # Compare diff signal in promoters and nearby enhancers
#' plot1 <- plotDiffScatter(object = ce_examples,
#'                contrasts =  c('contrast01', 'contrast01'),
#'                experimentNames = c('Peaks', 'Peaks'),
#'                filterAnnoCols = c('simple_annotation', 'simple_annotation'),
#'                filterValues = c('Promoter', 'Distal Intergenic'),
#'                assignFunc = GenomicRanges::nearest)
#'
#' # Compare diff signal in all peaks for two different contrasts
#' plot2 <- plotDiffScatter(object = ce_examples,
#'                contrasts =  c('contrast01', 'contrast02'),
#'                experimentNames = c('Peaks', 'Peaks'))
#'
plotDiffScatter <- function(object, experimentNames, contrasts,
                            filterAnnoCols = NULL, filterValues = NULL,
                            assignFunc = findOverlaps,
                            geneAnnoCols="nearest_TSS",
                            colorByFDR="mean") {
    if(sum(experimentNames %in% names(object)) != 2) {
        stop(sprintf('No experiment called %s in ChrawExperiment object. ',
                     experimentNames[!experimentNames %in% names(object)]))
    }
    if(!colorByFDR %in% c("exp1", "exp2", "mean", "none"))
        stop("The parameter `colorByFDR=` must be either 'exp1', 'exp2', 'mean' or 'none'")
    rnaseqFlagCheck <- list(
        metadata(experiments(object)[[experimentNames[1]]])$isrnaseq,
        metadata(experiments(object)[[experimentNames[2]]])$isrnaseq )
    rnaseqFlagCheck <- !vapply(rnaseqFlagCheck, is.null, logical(1))
    if(any(rnaseqFlagCheck)){
        scatter.df <- assignPeaks2GeneDE(
            object,
            experimentNames = c(experimentNames[!rnaseqFlagCheck], experimentNames[rnaseqFlagCheck]),
            geneAnnoCols = geneAnnoCols,
            filterAnnoCols = filterAnnoCols,
            filterValues = filterValues )
        scatter.df <- as.data.frame(scatter.df)
        contrasts <- c(contrasts[!rnaseqFlagCheck], contrasts[rnaseqFlagCheck])
        if( which(rnaseqFlagCheck) == 1 & colorByFDR == "exp1" )
            colorByFDR <- "exp2"
        else if( which(rnaseqFlagCheck) == 1 & colorByFDR == "exp2" )
            colorByFDR <- "exp1"
        if( !is.null(filterValues) ){
            filterValues <- c(filterValues, "")
        }
    }else{
        ## Check if any of the var vectors is too long
        lapply(c('experimentNames', 'contrasts', 'filterAnnoCols', 'filterValues'), function(varName) {
            if(length(eval(parse(text = varName))) > 2)
                stop(sprintf('Too many arguments given to %s. Please provide only two %s.',
                             varName, varName))
        })
        ## Combine rowData from diff/same experiment and filter for certain annotation
        scatter.df <- assignPeaks2Peaks(
            object = object,
            experimentNames = experimentNames,
            filterAnnoCols = filterAnnoCols,
            filterValues = filterValues,
            assignFunc = assignFunc)
        scatter.df <- as.data.frame(scatter.df)
    }
    ## Prepare padj column names
    padj.cn <- sprintf('%s_padj', contrasts)
    if(experimentNames[1] == experimentNames[2]) padj.cn[2] <- paste0(padj.cn[2], '.1')
    ## Check if diffAnalysis was performed for the defined contrasts
    if(sum(padj.cn %in% colnames(scatter.df)) != 2) {
        stop(sprintf('Contrast %s is not available for experiment %s.
                 Please run the testForDiffSignal function.',
                 contrasts[!padj.cn %in% colnames(scatter.df)],
                 experimentNames[!padj.cn %in% colnames(scatter.df)]))
    }
    for( i in seq_len(length(padj.cn)) ){
        scatter.df[,padj.cn[i]] <- ifelse( is.na(scatter.df[,padj.cn[i]]), 1, scatter.df[,padj.cn[i]] )
        scatter.df[,padj.cn[i]] <- pmax( scatter.df[,padj.cn[i]], 10^-40 )
    }
    padj <- "padj"
    if( colorByFDR == "mean" ){
        ## Compute average -log10(padj)
        scatter.df$padj <- (-log10(scatter.df[,padj.cn[1]]) +
                            -log10(scatter.df[,padj.cn[2]]))/2
    }else if( colorByFDR == "exp1" ){
        scatter.df$padj <- -log10(scatter.df[,padj.cn[1]])
    }else if( colorByFDR == "exp2" ){
        scatter.df$padj <- -log10(scatter.df[,padj.cn[2]])
    }else{
        padj <- NULL
    }
    ## Prepare xlab/ylab
    gglabs <- lapply(seq_len(length(contrasts)), function(ii) {
        if(!is.null(filterValues[ii])) {
            sprintf('%s %s log2(FC)', filterValues[ii], contrasts[ii])
        } else {
            sprintf('%s log2(FC)', contrasts[ii])
        }
    })
    ## PLOT: Scatterplot - contrast1 log2(FC) vs contrast2 log2(FC)
    log2fc.cn <- sprintf('%s_log2FoldChange', contrasts)
    if(experimentNames[1] == experimentNames[2]) log2fc.cn[2] <- paste0(log2fc.cn[2], '.1')
    ggplot(data = scatter.df,
           aes(x = !!sym(log2fc.cn[1]), y = !!sym(log2fc.cn[2]), col = padj)) +
        geom_point() + xlab(gglabs[1]) + ylab(gglabs[2]) +
        geom_hline(yintercept = 0, linetype = 'dashed') +
        geom_vline(xintercept = 0, linetype = 'dashed') +
        scale_color_gradient2(low = 'grey', mid = 'white', high = 'darkred',
                              midpoint = -log10(0.01), name = 'Mean\n-log10(FDR)') +
        ##scale_size(trans = 'reverse', name = 'log10(Dist)') +
        theme_cowplot()
}


