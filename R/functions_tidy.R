#' Function to tidy a given experiment from a Chraw MultiAssayExperiment.
#'
#' @param object A ChrawExperiment object
#' @param experimentName The name of the count experiment, normally added using the `addCountExperiment` function.
#' @param assayName Name of the assay (default: counts) which should be converted from wide to long format.
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
  assay_data <- assay(object[[experimentName]], assayName = assayName)
  assay_data <- as.data.frame(assay_data)
  assay_data$RowID <- rownames(object[[experimentName]])
  row_names <- assay_data$RowID

  sample_data <- as.data.frame(colData(object))

  assay_data_long <- reshape(assay_data, direction = "long", varying = names(assay_data),
                             v.names = "value", timevar = "sampleName", times = names(assay_data),
                             idvar = row_names, new.row.names = 1:(nrow(assay_data) * ncol(assay_data)),
                             sep = "_")
  assay_data_long$RowID <- row_names
  combined_data <- merge(assay_data_long, sample_data, by.x = "sampleName", by.y = "row.names", all.x = TRUE)
  combined_data[, "sampleName.y"] <- NULL # remove duplicate column

  return(combined_data)
}
