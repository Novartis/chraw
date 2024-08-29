#' Function to test enrichment of peak annotation in up/down regulated peaks.
#'
#' @description The function computes the enrichment of peak annotation in a
#' subset of differential peaks using a Fisher's exact test. For example,
#' the enrichment of chromHMM states in differential peaks compared with all peaks.
#'
#' @param object A ChrawExperiment object.
#' @param experimentName The name of the count experiment to test, normally added using the `addCountExperiment` function.
#' @param contrastName The name of the differential peak contrast used in the `testForDiffSignal` function.
#' @param annoName The name of the peak annotation added with `annotateExperimentRegions` function. This name will be searched in the colnames of the `rowData()` of the 'experimentName'.
#' @param padjThres Threshold for the adjusted pvalue from the diff peak analysis with `testForDiffSignal`. Default: 0.05.
#' @param log2FcThrs Threshold for the log2FC from the diff peak analysis. Default: 0.
#' @param foreground Define which peak set showed be tested, either "up" or "down" regulated peaks. Default: 'up'.
#' @param background Define which peak set should be used as background in the fisher's exact test, either "all" non diff peaks, "up" and "down". Default: 'all'.
#' @param ... Additional parameters passed to `fisher.test()`.
#'
#' @return A data.frame containing following columns:
#' \itemize{
#' \item{'annoName': Annotation name from the peak rowData.}
#' \item{'n_anno_fg_peaks': Number of foreground peaks overlapping the annotation.}
#' \item{'n_anno_bg_peaks': Number of background peaks overlapping the annotation.}
#' \item{'conf1': Lower boundry of the log2 odds confidence interval from `fisher.test`.}
#' \item{'conf2': Upper boundry of the log2 odds confidence interval from `fisher.test`.}
#' \item{'odds': Log2 odds ratio from the test statistic of `fisher.test`.}
#' \item{'pvalue': pvalue from `fisher.test`.}
#' }
#'
#' @seealso fisher.test
#'
#' @examples
#' data(ce_examples)
#' ce_examples <- rewrite_paths(ce_examples)
#'
#' ce_examples <- testForDiffSignal(
#'    object=ce_examples, experimentName ="Peaks",
#'    design = ~condition,
#'    contrasts = list(
#'    agonist_3h=c("condition", "agonist_3h", "CTRL_3h")))
#' ce_examples <- annotateExperimentRegions( ce_examples, "Peaks" )
#' ftResults <- enrichAnno(ce_examples,"Peaks",contrastName="agonist_3h",annoName="simple_annotation" )
#'
#' @export
#'
enrichAnno <- function( object, experimentName, contrastName, annoName, padjThres = 0.05,
                        log2FcThrs = 0, foreground = 'up', background = 'all', ... ) {

  # Check if comparison make sense
  if(foreground == background) stop(paste0("Foreground condition can't be the
                                           same as the background."))

  # Check if experimentName exists
  if(!(experimentName %in% names(object))) stop(paste0(experimentName, " was not found",
                                                    " in the ChrawExperiment."))

  # Parse enriched peaks from ChrawExperiment
  peaks <- rowData(object[[experimentName]])

  # Check if anno and contrast columns exist
  if( !sprintf('%s_padj', contrastName) %in% colnames(peaks) ) {
    stop(paste0("Diff peak analysis for ", contrastName, " is not existing. ",
                "Please run testForDiffSignal before."))
  }
  if( !annoName %in% colnames(peaks) ) {
    stop(paste0("Annotation column is not existing. ",
                "Please run annotateExperimentRegions before."))
  }

  # Compute annotation in foreground peaks
  fgSign <- ifelse(foreground == 'up', 1, -1)
  iFg <- which(peaks[,sprintf('%s_padj', contrastName)] < padjThres &
                 fgSign*peaks[,sprintf('%s_log2FoldChange', contrastName)] > log2FcThrs)
  nFgAnno <- table(peaks[iFg, annoName])

  # Compute annotation in background peaks
  if( background == 'all' ) {
    nBgAnno <- table(peaks[-iFg, annoName])
  } else {
    iBg <- which(peaks[,sprintf('%s_padj', contrastName)] < padjThres &
                   -fgSign*peaks[,sprintf('%s_log2FoldChange', contrastName)] > log2FcThrs)
    nBgAnno <-  table(peaks[iBg, annoName])
  }

  # Filter background peak overlap with foreground
  nBgAnno <- nBgAnno[names(nBgAnno) %in% names(nFgAnno)]

  # Perform fisher's exact test enrichment
  ftResults <- lapply(names(nBgAnno), function(an) {
    cont.table <- matrix(c(nFgAnno[an], nBgAnno[an],
                           sum(nFgAnno) - nFgAnno[an],
                           sum(nBgAnno) - nBgAnno[an]), nrow = 2)
    if(sum(is.na(cont.table))) cont.table[is.na(cont.table)] <- 0
    ftResult <- fisher.test(cont.table, ...)

    data.frame(annoName = an,
               n_anno_fg_peaks = nFgAnno[an],
               n_anno_bg_peaks = nBgAnno[an],
               conf1 = log2(ftResult$conf.int[1]),
               conf2 = log2(ftResult$conf.int[2]),
               odds = log2(ftResult$estimate),
               pvalue = ftResult$p.value)
  })

  ftResults <- do.call(rbind, ftResults)
  return(ftResults)
}
