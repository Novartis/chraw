#' Function to tidy a given experiment from a Chraw MultiAssayExperiment.
#'
#' @param object A ChrawExperiment object
#' @param experimentName The name of the count experiment, normally added
#'   using the `addCountExperiment` function.
#' @param assayName Name of the assay (default: counts) which should be
#'   converted from wide to long format.
#'
#' @return Tidy dataframe of experiment in ChrawExperiment object
#' @export
#'
#' @examples
#' data('ce_examples')
#'
#' ce_examples <- rewrite_paths(ce_examples)
#'
#' tidyChrawExperiment(ce_examples, experimentName = 'Peaks')
#' tidyChrawExperiment(ce_examples, experimentName = 'Peaks', assayName = 'cpms')
#'
tidyChrawExperiment <- function(object, experimentName, assayName = 'counts') {
  # Tidy given experiment and join with sample meta.data
  se <- object[[experimentName]]
  ## `assay(se, assayName = ...)` silently falls through to the first assay,
  ## so the assay has to be selected by name explicitly.
  if( !assayName %in% assayNames(se) ){
    stop(sprintf("The assay '%s' was not found in the experiment '%s'.",
                 assayName, experimentName))
  }
  assay_data <- as.data.frame(assays(se)[[assayName]])
  ## Only the sample columns are reshaped; RowID identifies the region and
  ## must not be treated as if it were another sample.
  sampleCols <- colnames(assay_data)
  assay_data$RowID <- rownames(se)

  sample_data <- as.data.frame(colData(object))

  assay_data_long <- reshape(assay_data, direction = "long",
                             varying = sampleCols,
                             v.names = "value", timevar = "sampleName",
                             times = sampleCols,
                             idvar = "RowID",
                             new.row.names = seq_len(nrow(assay_data) *
                                                     length(sampleCols)),
                             sep = "_")
  combined_data <- merge(assay_data_long, sample_data, by.x = "sampleName", by.y = "row.names", all.x = TRUE)
  combined_data[, "sampleName.y"] <- NULL # remove duplicate column

  return(combined_data)
}
