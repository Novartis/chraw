#' Functions to plot QC statistics
#' @rdname plottingFunctions
#' @description These functions parse the QC statistics from the json QC
#' file produced by the cromwell pipeline. These functions allow users to
#' plots the number of mapped and unmapped reads for each sample, number of
#' PCR duplicates per sample, fraction of reads in known annotation,
#' fraction of reads mapping to chrM, fraction of reads in peaks and
#' fragment length distributions.
#'
#' @param object A ChrawExperiment object.
#'
#' @examples
#' data(ce_examples)
#' ce_examples <- rewrite_paths(ce_examples)
#'
#' plot_ms <- plotMappingStats( ce_examples )
#' plot_pcr <- plotPCRDuplicateStats( ce_examples )
#'
#' @return A `ggplot` object.
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#' @export
plotMappingStats <- function( object ){
    if( !is( object, "ChrawExperiment" ) ){
        stop("The parameter 'object' must be a ChrawExperiment object")
    }
    qcData <- readQcJsonFromCE( object )
    mappingData <- do.call(rbind, lapply( names(qcData), function(x){
        do.call(rbind, lapply(names( qcData[[x]][["align"]][["samstat"]] ), function(y){
            nm <- paste(x, y, sep="_")
            dfs <- as.data.frame(qcData[[x]][["align"]][["samstat"]][[y]])
            dfs$samp <- nm
            dfs
        } ) )
    } ) )

    mappingData$unmapped_reads <- mappingData$total_reads - mappingData$mapped_reads
    mappingDataforPlot <- mappingData
    mappingData <- mappingData[,c("samp", "mapped_reads", "unmapped_reads")]
    mappingData <- reshape(mappingData, direction = "long",
                            varying = c("mapped_reads", "unmapped_reads"),
                            v.names = "read_number",
                            timevar = "type",
                            times = c("mapped_reads", "unmapped_reads"),
                            idvar = "samp",
                            sep = "_")

    # Remove any rows with missing values
    mappingData <- mappingData[complete.cases(mappingData), ]
    # Reset the row names of the reshaped dataframe
    row.names(mappingData) <- NULL
    mappingData$type <- gsub("_reads", "", mappingData$type )
    mappingData$type <- factor(
        ifelse( mappingData$type == "mapped", "Mapped", "Unmapped" ),
        levels=c("Unmapped", "Mapped") )

    mappingData %>%
        ggplot( aes( samp, read_number/2/1000000, fill=type ) ) +
        geom_bar(stat="identity") +
        geom_text(data = mappingDataforPlot, aes(x = samp, y = 3 + mapped_reads/2/10^6,
                                          label = paste0(pct_mapped_reads, '%')),
                  inherit.aes = FALSE) +
        themeChraw() +
        labs(y=chrawLabels$fragmentsM,
             x=chrawLabels$sample,
             fill="Alignment status",
             title="Read mapping",
             caption=paste("Read counts halved to fragments (paired-end).",
                           "Labels show the percentage of mapped reads."))
}

#' @rdname plottingFunctions
#' @export
plotPCRDuplicateStats <- function( object ){
    if( !is( object, "ChrawExperiment" ) ){
        stop("The parameter 'object' must be a ChrawExperiment object")
    }
    qcData <- readQcJsonFromCE( object )
    pcrDupData <- do.call(rbind, lapply( names(qcData), function(x){
        do.call(rbind, lapply(names(qcData[[x]][["align"]][["dup"]]), function(l){
            rs <- as.data.frame(qcData[[x]][["align"]][["dup"]][[l]])
            rs$samp <- paste(x, l, sep="_")
            rs
        }))
    }))
    pcrDupData$pct_duplicate_reads <- round(pcrDupData$pct_duplicate_reads, 2)
    if(sum(pcrDupData$paired_reads) == 0) {
      pcrDupData <- pcrDupData[, c("samp", "unpaired_reads", "unpaired_duplicate_reads", "pct_duplicate_reads")]

      # Create a new column "unpaired_reads" by subtracting unpaired_duplicate_reads from unpaired_reads
      pcrDupData$unpaired_reads <- pcrDupData$unpaired_reads - pcrDupData$unpaired_duplicate_reads

      # Reshape the data from wide to long format
      pcrDupData <- reshape(pcrDupData, direction = "long",
                              varying = c("unpaired_reads", "unpaired_duplicate_reads"),
                              v.names = "read_number",
                              timevar = "type",
                              times = c("unpaired_reads", "unpaired_duplicate_reads"),
                              idvar = "samp",
                              sep = "_")

      # Create a new column "type" based on the condition using grepl()
      pcrDupData$type <- ifelse(grepl("duplicate", pcrDupData$type), "duplicate", "singlicate")

      # Remove any rows with missing values
      pcrDupData <- pcrDupData[complete.cases(pcrDupData), ]

      # Reset the row names of the reshaped dataframe
      row.names(pcrDupData) <- NULL
      } else {
        pcrDupData <- pcrDupData[, c("samp", "paired_reads", "paired_duplicate_reads", "pct_duplicate_reads")]

        # Create a new column "paired_reads" by subtracting paired_duplicate_reads from paired_reads
        pcrDupData$paired_reads <- pcrDupData$paired_reads - pcrDupData$paired_duplicate_reads

        # Reshape the data from wide to long format
        pcrDupData <- reshape(pcrDupData, direction = "long",
                              varying = c("paired_reads", "paired_duplicate_reads"),
                              v.names = "read_number",
                              timevar = "type",
                              times = c("paired_reads", "paired_duplicate_reads"),
                              idvar = "samp",
                              sep = "_")

        # Create a new column "type" based on the condition using grepl()
        pcrDupData$type <- ifelse(grepl("duplicate", pcrDupData$type), "duplicate", "singlicate")

        # Remove any rows with missing values
        pcrDupData <- pcrDupData[complete.cases(pcrDupData), ]

        # Reset the row names of the reshaped dataframe
        row.names(pcrDupData) <- NULL
    }
    pcrDupData$type <- factor(
        ifelse( pcrDupData$type == "duplicate", "Duplicate", "Unique" ),
        levels=c("Duplicate", "Unique") )

    ## One label per sample, sitting just above that sample's own bar. Using
    ## the long table directly would place every label at the same height.
    labelData <- aggregate( read_number ~ samp, data=pcrDupData, FUN=sum )
    labelData$pct_duplicate_reads <- pcrDupData$pct_duplicate_reads[
        match( labelData$samp, pcrDupData$samp )]
    labelData$labelHeight <- labelData$read_number/10^6 +
        0.05*max(labelData$read_number)/10^6

    pcrDupData %>%
        ggplot( aes( samp, read_number/1000000, fill=type ) ) +
        geom_bar(stat="identity") +
        geom_text(data = labelData,
                  aes(x = samp,
                      y = labelHeight,
                      label = paste0(pct_duplicate_reads, '%')),
                  inherit.aes = FALSE) +
        themeChraw() +
        labs(y=chrawLabels$fragmentsM,
             x=chrawLabels$sample,
             fill="Duplicate status",
             title="PCR duplication",
             caption="Labels show the percentage of duplicated fragments.")
}

#' @rdname plottingFunctions
#' @export
plotChrMStats <- function( object ){
    if( !is( object, "ChrawExperiment" ) ){
        stop("The parameter 'object' must be a ChrawExperiment object")
    }
    qcData <- readQcJsonFromCE( object )
    chrMData <- do.call(rbind, lapply( names(qcData), function(x){
        do.call(rbind, lapply(names(qcData[[x]][["align"]][["frac_mito"]]), function(l){
            rs <- as.data.frame(qcData[[x]][["align"]][["frac_mito"]][[l]])
            rs$samp <- paste(x, l, sep="_")
            rs
        }) )
    } ) )
    if( is.null(chrMData) || nrow(chrMData) == 0 ){
        stop(paste("The QC report does not contain the mitochondrial read",
                   "statistics ('frac_mito'). These are produced by the",
                   "ATAC-seq pipeline only."))
    }
    chrMData$frac_mito_reads <- round(100*chrMData$frac_mito_reads, 2)
    chrMDataforPlot <- chrMData

    # Reshape the data from wide to long format
    chrMData <- reshape(chrMData, direction = "long",
                            varying = c("non_mito_reads", "mito_reads"),
                            v.names = "read_number",
                            timevar = "type",
                            times = c("non_mito_reads", "mito_reads"),
                            idvar = "samp",
                            sep = "_")

    # Create a new column "type" based on the condition using grepl()
    chrMData$type <- factor(
        ifelse(grepl("non_mito_reads", chrMData$type), "Other chromosomes", "chrM"),
        levels=c("chrM", "Other chromosomes") )

    # Reset the row names of the reshaped dataframe
    row.names(chrMData) <- NULL

    chrMData %>%
        ggplot( aes( samp, read_number/1000000, fill=type ) ) +
        geom_bar(stat="identity") +
        geom_text(data = chrMDataforPlot,
                  aes(x = samp,
                      y = 3 +  non_mito_reads/10^6,
                      label = paste0(frac_mito_reads, '%')),
                  inherit.aes = FALSE) +
        themeChraw() +
        labs(y=chrawLabels$readsM,
             x=chrawLabels$sample,
             fill="Chromosome",
             title="Mitochondrial read content",
             caption=paste("Counts are reads, not fragments.",
                           "Labels show the percentage of chrM reads."))
}

#' @rdname plottingFunctions
#' @export
plotFripStats <- function( object ){
    if( !is( object, "ChrawExperiment" ) ){
        stop("The parameter 'object' must be a ChrawExperiment object")
    }
    qcData <- readQcJsonFromCE( object )
    fripData <- lapply( names(qcData), function(x){
        lapply(names(qcData[[x]][["peak_enrich"]][["frac_reads_in_peaks"]])[1], function(l){
            ii <- grep(names(qcData[[x]][["peak_enrich"]][["frac_reads_in_peaks"]][[l]]),
                       pattern = '(rep[0-9]|pooled)$')
            rs <- as.data.frame(qcData[[x]][["peak_enrich"]][["frac_reads_in_peaks"]][[l]][ii])
            colnames(rs) <- names(qcData[[x]][["peak_enrich"]][["frac_reads_in_peaks"]][[l]])[ii]
            rs$samp <- paste(x, l, sep="_")
            rs
        })
    })
    ## As above: flatten the nested list before binding the data frames.
    fripData <- do.call(rbind, unlist(fripData, recursive=FALSE))

    if(is.null(fripData) || nrow(fripData) == 0) {
      stop('Provided QC JSON is not containing the FRIP statistics.')
    }

    ## Carry the replicate names over into the reshaped table, so that the
    ## sample labels read "<sample>_<peakCaller>_rep1" and not "..._1".
    repNames <- colnames(fripData)[seq_len(ncol(fripData) - 1)]
    fripData <- reshape(fripData, direction = "long",
                              varying = seq_len(ncol(fripData) - 1),
                              v.names = "frip",
                              timevar = "rep",
                              times = repNames,
                              idvar = "samp",
                              new.row.names = NULL,
                              sep = "_")

    # Combine column values of "samp" and "rep" using paste()
    fripData$samp <- paste(fripData$samp, fripData$rep, sep = "_")

    # Remove rows with missing values in "frip"
    fripData <- subset(fripData, !is.na(frip))
    fripData %>%
        ggplot( aes( samp, 100*frip ) ) + geom_bar(stat="identity") +
        themeChraw() +
        labs(y="Reads in peaks (%)",
             x=chrawLabels$sample,
             title="Fraction of reads in peaks (FRiP)")
}

#' @rdname plottingFunctions
#' @export
plotFracReadsInAnnot<- function( object ){
    if( !is( object, "ChrawExperiment" ) ){
        stop("The parameter 'object' must be a ChrawExperiment object")
    }
    qcData <- readQcJsonFromCE( object )
    annoData <- lapply( names(qcData), function(x){
        lapply(names(qcData[[x]][["align"]][["frac_reads_in_annot"]]), function(l){
            rs <- as.data.frame(qcData[[x]][["align"]][["frac_reads_in_annot"]][[l]])
            colnames(rs) <- names(qcData[[x]][["align"]][["frac_reads_in_annot"]][[l]])
            rs$samp <- paste(x, l, sep="_")
            rs
        })
    })
    ## The nested lapply above yields a list of lists; flatten it first,
    ## otherwise rbind() stacks the lists into a matrix instead of binding
    ## the data frames.
    annoData <- do.call(rbind, unlist(annoData, recursive=FALSE))
    if( is.null(annoData) || nrow(annoData) == 0 ){
        stop(paste("The QC report does not contain the fraction of reads in",
                   "annotated regions ('frac_reads_in_annot'). These are",
                   "produced by the ATAC-seq pipeline only."))
    }
    ## Keep the annotation names ("fri_dhs", "fri_prom", ...) instead of the
    ## column indices, otherwise the legend keys come out as 1, 2, 3, ...
    annoNames <- colnames(annoData)[seq_len(ncol(annoData) - 1)]
    annoData <- reshape(annoData, direction = "long",
                            varying = seq_len(ncol(annoData) - 1),
                            v.names = "freq",
                            timevar = "anno",
                            times = annoNames,
                            idvar = "samp",
                            new.row.names = NULL,
                            sep = "_")

    # Remove 'fri_' from the "anno" column using gsub()
    annoData$anno <- gsub("fri_", "", annoData$anno)
    ## Spell out the ENCODE abbreviations for the legend.
    prettyAnno <- c( dhs="DNase hypersensitive sites", prom="Promoters",
                     enh="Enhancers", blacklist="ENCODE blacklist" )
    annoData$anno <- ifelse( annoData$anno %in% names(prettyAnno),
                             prettyAnno[annoData$anno], annoData$anno )
    annoData %>%
        ggplot( aes( samp, 100*freq, fill = anno) ) +
        geom_bar(stat="identity", position = 'dodge') +
        themeChraw() +
        labs(y="Reads in annotation (%)",
             x=chrawLabels$sample,
             fill="Genomic annotation",
             title="Reads overlapping annotated regions",
             caption=paste("Categories overlap, so the bars of a sample do",
                           "not add up to 100%."))
}
#' @rdname plottingFunctions
#' @export
plotFragmentLengthDist <- function( object ){
  if(!is( object, "ChrawExperiment" )){
    stop("The parameter 'object=' must be a ChrawExperiment object")
  }
  if( !"FragLengthDist" %in% names( experiments( object ) ) ){
    stop("Fragment length distribution data not found. Please run `addFragmentLengthDist()` before running `plotFragLengthDist()`")
  }
  if( !"sampleName" %in% colnames(colData(object))){
    stop("'sampleName' should be the name of the samples column. Change the name of this column to run this function")
  }
  plottingData <- as.data.frame(experiments(object)[["FragLengthDist"]])
  plottingData$fragmentSize <- rownames(plottingData)
  numTimes <- ncol(plottingData) - 1

  # Reshape the data from wide to long format
  plottingData <- reshape(plottingData, direction = "long",
                          varying = names(plottingData)[seq_len(numTimes)],
                          v.names = "frequency",
                          times = names(plottingData)[seq_len(numTimes)],
                          idvar = "fragmentSize",
                          new.row.names = NULL)
  # Convert the necessary columns from `colData(object)` to a data frame
  colData <- data.frame(sampleName = colData(object)$sampleName,
                        cond = colData(object)$condition,
                        check.names = FALSE)
  # Merge the `plottingData` dataframe with the `colData` dataframe
  plottingData <- merge(plottingData, colData, by.x = "time", by.y = "sampleName", all.x = TRUE)
  # Rename the columns "sampleName" and "cond" to "sampleID" and "sampleGroup"
  plottingData <- transform(plottingData, sampleID = time, sampleGroup = cond)
  plottingData$time <- NULL
  plottingData$cond <- NULL

  plottingData$rep <- gsub( "\\S+_(\\S+)$", "rep\\1", plottingData$sampleID ) ## strong assumption on name structures
  plottingData$fragmentSize <- as.numeric( plottingData$fragmentSize )
  plottingData %>%
    ggplot( aes( fragmentSize, frequency, color=sampleGroup, group=sampleID, linetype=rep ) ) +
    stat_smooth(geom="line", method="loess", span=0.02, se=FALSE, alpha=0.7, linewidth=0.5) +
    facet_wrap( ~sampleGroup ) +
    themeChraw( rotateXLabels=FALSE ) +
    ## The facet strip already names the condition, so only the replicate
    ## legend is worth keeping: without it the line types cannot be read.
    guides( color="none" ) +
    labs(x=chrawLabels$fragLength,
         y="Fraction of fragments",
         linetype="Replicate",
         title="Fragment length distribution",
         caption="Curves are LOESS-smoothed (span = 0.02).")
}

plotPCA.ChrawExperiment <- function( object, experimentName, prinComps=c(1, 2), colourGrouping, includeSamples=NULL, transformation=log10, ntop=2000  ){
    if( !experimentName %in% names( experiments(object) ) ){
        stop(sprintf("The experiment '%s' was not found in the ChawExperiment object", experimentName))
    }
    if( !is( experiments(object)[[experimentName]], "SummarizedExperiment" ) ){
        stop(sprintf("The experiment must be a SummarizedExperiment added with the `addCountExperiment()` function") )
    }
    se <- experiments(object)[[experimentName]]
    if(!is.null(includeSamples)){
        if( !all( includeSamples %in% colnames( se ) ) )
            stop(sprintf("The sample(s) %s was(were) not found in the ChrawExperiment object",
                         paste(includeSamples[!includeSamples %in% colnames( se )], collapse=", ") ) )
        se <- se[,includeSamples]
    }
    if( !all(prinComps >= 1 & prinComps <= ncol(se))){
        stop( sprintf("The specified PC loadings are out of range, please specify values in between %d and %d",
                     1, ncol(se) ) )
    }
    if( !colourGrouping %in% colnames( colData(object) ) ){
        stop(sprintf("The variable '%s' not found in the ChrawExperiment object. The variable must be a column of `colData(object)`.",
                     colourGrouping))
    }
    countData <- assays(se)[["counts"]]
    sizeFacs <- estimateSizeFactorsForMatrix( countData )
    countData <- sweep(countData, 2, sizeFacs, FUN=`/`)
    countData <- transformation( countData )
    rv <- rowVars( countData )
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca <- prcomp( t( countData[select, ] ) )
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    pcaDF <- pca$x[,prinComps]
    pcaDF <- cbind( pcaDF,
                   as.data.frame(colData(object)[rownames(pcaDF),,drop=FALSE] ))
    xvar <- paste0("PC", prinComps[1])
    yvar <- paste0("PC", prinComps[2])
    pcaDF %>%
        ggplot( aes( get(xvar), get(yvar), col=get(colourGrouping) ) ) +
        geom_point() +
        xlab(sprintf("%s (%s%% of variance)", xvar,
                     round(percentVar[prinComps[1]] * 100))) +
        ylab(sprintf("%s (%s%% of variance)", yvar,
                     round(percentVar[prinComps[2]] * 100))) +
        themeChraw( rotateXLabels=FALSE ) +
        labs(color=colourGrouping,
             title=sprintf("PCA of experiment '%s'", experimentName),
             subtitle=sprintf(
                 "%d most variable regions, size-factor normalised counts",
                 min(ntop, nrow(countData)))) +
        coord_fixed()
}

#' Plots to explore the sample-to-sample correlations between samples.
#'
#' These functions document the S4 method plotPCA for
#' count experiments of a ChrawExperiment object
#'
#' @docType methods
#' @name plotPCA
#' @rdname plottingSampCors
#' @aliases plotPCA plotPCA,ChrawExperiment-method
#'
#' @param object A ChrawExperiment object.
#' @param experimentName The name of the count experiment, normally added
#'   using the `addCountExperiment` function.
#' @param prinComps A numeric vector of length 2 that specifies the PC
#'   loadings to plot.
#' @param colourGrouping A variable name from the metadata to color the
#'   points. Must be a column of `colData(object)`.
#' @param includeSamples A character vector specifying which samples to
#'   include, these must be in `rownames(colData(object))`.
#' @param transformation A function specifying a function to apply to the
#'   data. We recommend a variance-stabilizing transformation.
#' @param ntop A numeric value indicating the number of top features to use
#'   for principal components, selected by highest variance.
#'
#' @importFrom DESeq2 estimateSizeFactorsForMatrix
#' @importFrom matrixStats rowVars
#' @importMethodsFrom BiocGenerics plotPCA
#'
#' @return A `ggplot` object with the samples projected on the requested
#' principal components.
#'
#' @examples
#' data(ce_examples)
#'
#' ce_examples <- rewrite_paths(ce_examples)
#'
#' plot <- plotPCA( ce_examples, "Peaks", colourGrouping="condition", transformation=asinh )
#'
#' @export
setMethod( "plotPCA", signature(object="ChrawExperiment"), plotPCA.ChrawExperiment )


#' Plots to explore the sample-to-sample correlations between samples.
#'
#' This function plots the pairwise correlation between samples for
#' count experiment contained within a ChrawExperiment object.
#'
#' @name plotSampleCorrelations
#'
#' @param object A ChrawExperiment object.
#' @param experimentName The name of the count experiment, normally added
#'   using `addCountExperiment` function.
#' @param annotationColumns A variable name from the metadata to color the
#'   points. Must be a column of colData(object).
#' @param includeSamples A character vector specifying which samples to
#'   include. These must be in `rownames(colData(object))`.
#' @param transformation A function specifying a function to apply to the
#'   data. We recommend a variance-stabilizing transformation.
#' @param ntop A numeric value indicating the number of top features to use
#'   for principal components, selected by highest variance.
#' @param annotationColors A list of colors passed to the
#'   `HeatmapAnnotation()` function. For details check `?HeatmapAnnotation`.
#' @param sampleLabelColumn A column name of the `colData(object)` to use for
#'   labelling. If missing, the default sample names will be used.
#' @param ... Further parameters passed to the `Heatmap()` function. Check
#'   `?Heatmap` for further details.
#'
#' @importFrom circlize colorRamp2
#' @importFrom ComplexHeatmap HeatmapAnnotation Heatmap
#' @importMethodsFrom BiocGenerics plotPCA
#'
#' @return A `Heatmap` object from the ComplexHeatmap package, showing the
#' Spearman correlation between every pair of samples.
#'
#' @examples
#' data(ce_examples)
#' ce_examples <- rewrite_paths(ce_examples)
#'
#' plot <- plotSampleCorrelations( ce_examples, "Peaks", annotationColumns="condition" )
#'
#' @export
plotSampleCorrelations <- function( object, experimentName, annotationColumns=NULL,
                                   transformation=log10, includeSamples=NULL, ntop=2000, annotationColors=NULL,
                                   sampleLabelColumn=NULL, ... ){

    if( !experimentName %in% names( experiments(object) ) ){
        stop(sprintf("The experiment '%s' was not found in the ChawExperiment object", experimentName))
    }
    if( !is( experiments(object)[[experimentName]], "SummarizedExperiment" ) ){
        stop(sprintf("The experiment must be a SummarizedExperiment added with the `addCountExperiment()` function") )
    }
    se <- experiments(object)[[experimentName]]
    if(!is.null(includeSamples)){
        if( !all( includeSamples %in% colnames( se ) ) )
            stop(sprintf("The sample(s) %s was(were) not found in the ChrawExperiment object",
                         paste(includeSamples[!includeSamples %in% colnames( se )], collapse=", ") ) )
        se <- se[,includeSamples]
    }
    if( !all(annotationColumns %in% colnames( colData(object) )) ){
        stop(sprintf("The variables '%s' not found in the ChrawExperiment object. The variable must be a column of `colData(object)`.",
                     annotationColumns[!annotationColumns %in% colnames( colData(object) )] ) )
    }
    ## se <- experiments(object)[[experimentName]]
    countData <- assays(se)[["counts"]]
    oldNames <- colnames(countData)
    if( !is.null(sampleLabelColumn) ){
        if( !sampleLabelColumn %in% colnames( colData(object) )){
            stop(
                sprintf("The variables '%s' not found in the ChrawExperiment object. The `sampleLabelColumn=` parameter must be a column of `colData(object)`.",
                        sampleLabelColumn) )
        }
        newNames <- colData(object)[oldNames,sampleLabelColumn]
        colnames(countData) <- newNames
    }else{
        newNames <- oldNames
    }
    sizeFacs <- estimateSizeFactorsForMatrix( countData )
    countData <- sweep(countData, 2, sizeFacs, FUN=`/`)
    countData <- transformation( countData )
    rv <- rowVars( countData )
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    corMat <- cor( ( countData[select, ] ), method="spearman" )
    colFun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    if( !is.null(annotationColumns) ){
        annoDf <- as.data.frame(colData(object)[oldNames,annotationColumns,drop=FALSE])
        rownames(annoDf) <- newNames
        ## there must be a better solution to handle passing missing arguments ##
        if( !is.null(annotationColors) ){
            topAnn <- HeatmapAnnotation( df=annoDf, col=annotationColors )
        }else{
            topAnn <- HeatmapAnnotation( df=annoDf )
        }
    }else{
        topAnn <- NULL
    }
    Heatmap( corMat,
            name="Spearman\ncorrelation",
            col=colFun,
            top_annotation=topAnn,
            column_title=sprintf("Sample correlation of experiment '%s'",
                                 experimentName),
            row_title="Samples",
            column_names_rot=45,
            heatmap_legend_param=list(at=c(-1, -0.5, 0, 0.5, 1)),
            ... )
}


plotMA.ChrawExperiment <- function( object, experimentName, contrastName=NULL, FDR=0.1, meanFilter=0 ){
    diffRes <- pullDiffResults( object, experimentName=experimentName, contrastName=contrastName, outputFormat="df_long" )
    diffRes$padj <- ifelse(is.na(diffRes$padj), 1, diffRes$padj)
    diffRes$significant <- factor( ifelse( diffRes$padj < FDR, "yes", "no") )
    diffRes <- diffRes[diffRes$baseMean > meanFilter, ]
    diffRes %>%
        ggplot( aes( log10(baseMean), log2FoldChange, col=significant ) ) +
        geom_point() +
        facet_wrap( ~contrastName ) +
        geom_hline(yintercept=0, col="#4393c3", linewidth=1.2) +
        scale_color_manual(values=c( yes="#b2182b80", no="#4d4d4d80" ),
                           labels=c( yes="Yes", no="No" )) +
        themeChraw( rotateXLabels=FALSE ) +
        labs(x=expression("Mean normalised count (" * log[10] * ")"),
             y=chrawLabels$log2fc,
             col=sprintf("FDR < %s", FDR),
             title="MA plot",
             caption=sprintf(
                 "DESeq2 base mean and fold change; regions with a base mean below %s are not shown.",
                 meanFilter))
}

#' Plots to explore the sample-to-sample correlations between samples.
#'
#' These functions document the S4 method plotPCA for
#' count experiments of a ChrawExperiment object
#'
#' @docType methods
#' @name plotMA
#' @rdname plotDiff
#' @aliases plotMA plotMA,ChrawExperiment-method
#'
#' @param object A ChrawExperiment object.
#' @param experimentName The name of the count experiment, normally added
#'   using the `addCountExperiment` function.
#' @param contrastName A name of the contrast name to extract. For a list of
#'   the contrast names, consider running `getDiffResultSummary()`. If missing,
#'   all available results will be returned.
#' @param FDR An FDR threshold to assess significance.
#' @param meanFilter A numeric value to filter from plotting peaks below
#'   these threshold.
#'
#' @return A `ggplot` object, faceted by contrast.
#'
#' @examples
#' data(ce_examples)
#' ce_examples <- rewrite_paths(ce_examples)
#'
#' plotMA( ce_examples, experimentName = "Peaks" )
#' plotVolcano( ce_examples, experimentName = "Peaks" )
#'
#' @importFrom BiocGenerics plotMA
#' @export
setMethod( "plotMA", signature(object="ChrawExperiment"), plotMA.ChrawExperiment )

plotVolcano.ChrawExperiment <- function( object, experimentName, contrastName=NULL, FDR=0.1, meanFilter=0 ){
    diffRes <- pullDiffResults( object, experimentName=experimentName, contrastName=contrastName, outputFormat="df_long" )
    diffRes$padj <- ifelse(is.na(diffRes$padj), 1, diffRes$padj)
    diffRes$significant <- factor( ifelse( diffRes$padj < FDR, "yes", "no") )
    diffRes <- diffRes[diffRes$baseMean > meanFilter,]
    diffRes %>%
        ggplot( aes( log2FoldChange, -log10( pvalue ), col=significant ) ) +
        geom_point() +
        facet_wrap( ~contrastName ) +
        geom_vline(xintercept=0, col="#4393c3", linewidth=1.2) +
        scale_color_manual(values=c( yes="#b2182b80", no="#4d4d4d80" ),
                           labels=c( yes="Yes", no="No" )) +
        themeChraw( rotateXLabels=FALSE ) +
        labs(x=chrawLabels$log2fc,
             y=expression(-log[10] ~ "(raw p-value)"),
             col=sprintf("FDR < %s", FDR),
             title="Volcano plot",
             caption=paste("The y axis shows the raw p-value;",
                           "the colour marks the adjusted p-value (FDR)."))
}

#' @name plotVolcano
#' @rdname plotDiff
#' @aliases plotVolcano plotVolcano,ChrawExperiment-method
#' @export
setMethod( "plotVolcano", signature(object="ChrawExperiment"), plotVolcano.ChrawExperiment )

#' Functions for signal meta profiles.
#' @rdname metaProfiles
#' @description These functions provide functionalities to plot the
#' signal meta profiles for selected genomic regions. These functions
#' provide options to either return the meta profiles as R objects for further
#' analysis or plot them as heatmaps.
#'
#' @param object A ChrawExperiment object.
#' @param regions GRanges objects containing the regions to plot.
#' @param annotationColumns A variable name from the metadata. Must be a
#'   column of colData(object).
#' @param includeSamples A character vector specifying which samples to
#'   include, these must be in `rownames(colData(object))`.
#' @param inputColumn A column of `colData(object)` with the paths to the
#'   files used to calculate the metaprofiles. These could be either
#'   `bigwigFile` or `bamFile`.
#' @param outputFormat A character vector specifying the output format of the
#'   metaprofiles. These could be either "ScoreMatrixList" for a
#'   ScoreMatrixList object of the package genomation, "list" for a normal R
#'   list or "df_long" for a long-formatted data frame.
#' @param ... Parameters passed to the `ScoreMatrixList()` function.
#'
#' @return For `calculateMetaProfiles()`, a `ScoreMatrixList`, a list of
#' matrices or a long data frame, according to `outputFormat`. For
#' `reformatMatrixScores()`, a long data frame or a list of matrices.
#' For `plotMetaProfileHeatmaps()`, a `HeatmapList` object.
#'
#' @importFrom genomation ScoreMatrixList
#' @export
#'
#' @examples
#' data(ce_examples)
#'
#' ce_examples <- rewrite_paths(ce_examples)
#'
#' peaks <- GenomicRanges::resize(rowRanges(ce_examples[["Peaks"]])[seq_len(100),], 10^3, fix = 'center')
#' matrix <- calculateMetaProfiles(object = ce_examples, regions = peaks)
#'
calculateMetaProfiles <- function( object, regions, annotationColumns=NULL,
                                  includeSamples=NULL, inputColumn="bigwigFile",
                                  outputFormat="ScoreMatrixList",  ... ){
    stopifnot(is( object, "ChrawExperiment" ))
    cDat <- colData( object )
    if(!is.null(includeSamples)){
        if( !all( includeSamples %in% rownames(cDat))){
            stop(sprintf("The sample(s) %s was(were) not found in the ChrawExperiment object",
                         paste(includeSamples[!includeSamples %in% rownames(colData( object ))], collapse=", ") ) )
        }
        cDat <- cDat[rownames(cDat) %in% includeSamples,]
    }
    if(!outputFormat %in% c("ScoreMatrixList", "list", "long_df")){
        stop("The parameter 'outputFormat' must be either 'ScoreMatrixList', 'long_df' or 'list'")
    }
    if( !inputColumn %in% colnames(cDat) ){
        stop(sprintf("Invalid `inputColumn=` parameter. The '%s' variable must be a column of `colData(object)`.", inputColumn))
    }
    rs <- withBigwigPaths( cDat[,inputColumn], function( targets )
        ScoreMatrixList( targets=targets, windows=regions, ...  ) )
    names(rs) <- rownames(cDat)
    rs <- lapply( rs, function( x ){
        rownames(x) <- names(regions)
        x
    } )
    rs <- as(rs, "ScoreMatrixList")
    if( outputFormat != "ScoreMatrixList" ){
        rs <- reformatMatrixScores( rs, object, annotationColumns, outputFormat=outputFormat )
    }
    return(rs)
}

#' @rdname metaProfiles
#' @param matrixScores An object "ScoreMatrixList" obtained running the
#'   function `calculateMetaProfiles()`.
#' @export
reformatMatrixScores <- function( matrixScores, object, annotationColumns=NULL, outputFormat="long_df" ){
  stopifnot( is( matrixScores, "ScoreMatrixList" ) )
  stopifnot( is( object, "ChrawExperiment" ) )
  if(!outputFormat %in% c("list", "long_df")){
    stop("The parameter 'outputFormat' must be either 'long_df' or 'list'")
  }
  if( outputFormat == "list" ){
    rs <- lapply( matrixScores, function(x){ x@.Data } )
    names(rs) <- names(matrixScores)
  }else{
    rs <- lapply( names(matrixScores), function(x){
      #x <- names(matrixScores)[1]
      xname <- x
      x <- matrixScores[[x]]
      x <- x@.Data
      x <- as.data.frame( x )
      colnames(x) <- seq_len(ncol(x))
      x <- as.data.frame(x)
      x$regionID <- rownames(x)
      cols_to_exclude <- grep("^regionID$", colnames(x), value = TRUE, invert = TRUE)

      x <- reshape(x, direction="long",varying = cols_to_exclude, v.names = "score",
                   times = cols_to_exclude, timevar = "relative_pos", idvar = "regionID")
      #x <- x[order(x$regionID), ]
      rownames(x) <- NULL
      x$sampleName <- xname
      x$relative_pos <- as.numeric(as.character(x$relative_pos))
      x
    } )
    rs <- do.call(rbind, rs)
    if( !is.null(annotationColumns) ){
      if( !all(annotationColumns %in% colnames( colData(object) )) ){
        stop(sprintf("The variable(s) '%s' was(were) not found in the ChrawExperiment object. The values of `annotationColumns=` must be column names of `colData(object)`.",
                     annotationColumns[!annotationColumns %in% colnames( colData(object) )] ) )
      }
      cDat <- as.data.frame(colData(object)[,annotationColumns,drop=FALSE]) #%>%
      cDat$sampleName <- rownames(cDat)
      rs <- merge(rs, cDat, by = "sampleName", all.x = TRUE)
    }
  }
  return(rs)
}

#' @rdname metaProfiles
#' @param xcoords A vector of three strings with the labels of the left-most,
#'   centre and right-most positions of the heatmap. If `NULL` (the default),
#'   the labels are derived from the width of `regions` and are given in base
#'   pairs relative to the centre of the region.
#' @param sampleLabelColumn A column name of the `colData(object)` to use for
#'   labelling. If missing, the default sample names will be used.
#' @param keep Lower and upper quantile values to trim the matrix heatmap
#'   when plotting. Useful to avoid artifacts due to outliers.
#' @return A `Heatmap` object from the ComplexHeatmap package.
#'
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation anno_lines
#' @export
plotMetaProfileHeatmaps <- function( object, regions, xcoords=NULL,
                                    includeSamples, inputColumn="bigwigFile",
                                    sampleLabelColumn, keep=c(0, 0.99), ... ){
    if( !is.null(sampleLabelColumn) ){
        if( !sampleLabelColumn %in% colnames( colData(object) )){
            stop(
                sprintf("The variables '%s' not found in the ChrawExperiment object. The `sampleLabelColumn=` parameter must be a column of `colData(object)`.", sampleLabelColumn) )
        }
    }
    metaProfiles <-
        calculateMetaProfiles(
            object, regions,
            includeSamples=includeSamples,
            inputColumn=inputColumn,
            outputFormat="list", ... )
    stopifnot(names(metaProfiles) == rownames(colData(object)[includeSamples,]))
    if( !is.null(sampleLabelColumn) ){
        names(metaProfiles) <- colData(object)[includeSamples,sampleLabelColumn]
    }
    rMeans <- lapply( metaProfiles, function(x){
        rowMeans(x)
    } )
    rOrder <-  order(rowMeans(do.call(cbind, rMeans)), decreasing=TRUE)
    ## Label the x axis with the distance to the centre of the regions rather
    ## than with bare numbers, unless the user asked for their own labels.
    if( missing(xcoords) || is.null(xcoords) ){
        halfWidth <- round( median( width( regions ) ) / 2 )
        xcoords <- c( sprintf("-%s bp", halfWidth), "centre",
                      sprintf("+%s bp", halfWidth) )
    }
    ## A single colour scale across all samples, so that the panels can be
    ## compared with each other, and one shared legend that names the quantity.
    allValues <- unlist( metaProfiles, use.names=FALSE )
    globalBounds <- quantile( allValues, keep )
    col_fun <- colorRamp2( globalBounds, c("white", "red") )
    signalName <- ifelse( inputColumn == "bigwigFile",
                         "Coverage\n(bigwig signal)", "Coverage\n(reads)" )
    hList <- lapply( seq_along(metaProfiles), function(i){
        x <- names(metaProfiles)[i]
        mat <- metaProfiles[[x]]
        clMeans <- colMeans(mat)
        mat <- pmax( mat, globalBounds[1] )
        mat <- pmin( mat, globalBounds[2] )
        colnames(mat) <- rep("", ncol(mat))
        colnames(mat)[1L] <- xcoords[1]
        colnames(mat)[ncol(mat)] <- xcoords[3]
        colnames(mat)[round(ncol(mat)/2)] <- xcoords[2]
        ha <- HeatmapAnnotation(mean=anno_lines(clMeans, axis=TRUE),
                               annotation_label="Mean signal",
                               show_annotation_name=(i == 1L) )
        ## ComplexHeatmap requires a unique name per panel, so the shared
        ## colour scale is named through the legend title of the first panel
        ## instead of through `name`.
        Heatmap( mat[rOrder,], name=x, col = col_fun,
                cluster_columns = FALSE,
                cluster_rows=FALSE, show_row_names=FALSE,
                show_heatmap_legend = (i == 1L),
                heatmap_legend_param = list(title = signalName),
                row_title=sprintf("%d regions, sorted by mean signal",
                                  nrow(mat)),
                top_annotation=ha, column_title=x  )
    } )
    Reduce(`+`, hList)
}

#' @rdname enrichAnno
#' @param enrichRes A data.frame contating the annotation enrichment from
#'   `enrichAnno`.
#' @param pvalueThres P-value threshold to adjust the color scale.
#' @export
plotEnrichResults <- function( enrichRes, pvalueThres = 0.05 ) {
    if( sum(colnames(enrichRes) %in% c('annoName',  'conf1', 'conf2', 'odds', 'pvalue')) != 5)
        stop(paste0('Input data.frame is missing required columns. ',
                    'Please generate a new one with enrichAnno().'))
    ggplot(data = enrichRes, aes(x = reorder(annoName, -log10(pvalue)), y = odds,
                                 col = -log10(pvalue))) +
        geom_point() +
        geom_errorbar(aes(ymin = conf1, ymax = conf2), width = 0.15) +
        geom_hline(yintercept = 0, linetype = 'dashed') +
        scale_color_gradient2(low = 'grey', mid = 'white', high = 'darkred',
                              midpoint = -log10(pvalueThres),
                              name = expression(-log[10] ~ "(p-value)")) +
        labs(x = "Annotation category",
             y = expression(log[2] ~ "odds ratio (95% CI)"),
             title = "Annotation enrichment",
             caption = paste("Fisher's exact test; the colour shows the raw",
                             "p-value, not the FDR.")) +
        coord_flip() +
        themeChraw( rotateXLabels = FALSE )
}


#' Vplot functionality
#'
#' @description This function calculates the 'vplot' profile given a
#' ChrawExperiment object and a set of genomic ranges.
#'
#' @param object A ChrawExperiment object.
#' @param regions GRanges objects containing the regions to plot.
#' @param windowSize Window size in bp, used for resizing the regions.
#' @param returnMatrix Logical, indicating whether the vplot matrix should be
#'   returned instead of the plot.
#' @param sampleLabelColumn A character string with the column name of
#'   `colData(object)`. The values of this column will be used as sample
#'   labels.
#' @param includeSamples A character vector specifying which samples to
#'   include. These must be in `rownames(colData(object))`.
#' @param BPPARAM A BiocParallel instance.
#' @param zlim Numeric value between 0 and 1. Specifies the upper quantile
#'   value to be considered as maximum value in the color scale. Larger values
#'   in the heatmap are set to this max value (using pmin).
#' @param rowScaled Logical flag, indicating whether each rows of the
#'   heatmaps should be standarized (z-scores) or not.
#' @return Either a long formatted matrix, or a ggplot object.
#' @examples
#' data(ce_examples)
#' ce_examples <- rewrite_paths(ce_examples)
#'
#' indexBam( ce_examples )
#'
#' regions <- rowRanges(ce_examples[["Peaks"]])
#' vplot <- plotVProfile( ce_examples, regions = regions )
#'
#' @importFrom GenomicAlignments readGAlignmentPairs first second
#' @export
plotVProfile <- function( object, regions, windowSize=600, returnMatrix=FALSE, sampleLabelColumn = NULL, includeSamples=NULL, BPPARAM=SerialParam(), zlim=0.95, rowScaled=FALSE ){
  stopifnot(is( object, "ChrawExperiment" ))
  if( !is.null(sampleLabelColumn) &&
      !sampleLabelColumn %in% colnames( colData(object) )){
    stop(sprintf("The variables '%s' not found in the ChrawExperiment object. The `sampleLabelColumn=` parameter must be a column of `colData(object)`.", sampleLabelColumn))
  }
  possibleSamples <- rownames(colData(object))
  if(!is.null(includeSamples)){
    if( !all( includeSamples %in% possibleSamples ) )
      stop(sprintf("The sample(s) %s was(were) not found in the ChrawExperiment object",
                   paste(includeSamples[!includeSamples %in% possibleSamples], collapse=", ") ) )
  }else{
    includeSamples <- possibleSamples
  }
  regions <- resize(regions, windowSize, fix="center")
  flag <- scanBamFlag(isPaired = TRUE, isProperPair = TRUE,
                      isUnmappedQuery = FALSE, hasUnmappedMate = FALSE)
  param <- ScanBamParam(which = regions,
                        flag = flag)
  vMatrices <- bplapply( includeSamples, function(x){
    bf <- colData(object)[x,"bamFile"]
    alnp <- readGAlignmentPairs(bf, param = param, use.names = TRUE)
    leftMost <- ifelse( start(GenomicAlignments::first(alnp)) < end(GenomicAlignments::second(alnp)), start(GenomicAlignments::first(alnp)), end(GenomicAlignments::second(alnp)) )
    rightMost <- ifelse( start(GenomicAlignments::first(alnp)) < end(GenomicAlignments::second(alnp)), end(GenomicAlignments::second(alnp)), start(GenomicAlignments::first(alnp)) )
    readRanges <- GRanges(as.character(seqnames(GenomicAlignments::first(alnp))), IRanges(start=leftMost, end=rightMost))
    fragmentMidpoint <- resize(readRanges, 1, fix="center")
    mcols(fragmentMidpoint)$fragmentSize <- width(readRanges)
    ovls <- findOverlaps(fragmentMidpoint, regions)
    mcols(fragmentMidpoint)$relPos <- NA
    mcols(fragmentMidpoint)$relPos[queryHits(ovls)] <-
      start(fragmentMidpoint)[queryHits(ovls)] - start(regions)[subjectHits(ovls)] + 1
    fragmentMidpoint <- fragmentMidpoint[!is.na(mcols(fragmentMidpoint)$relPos)]
    if( length( fragmentMidpoint ) == 0 ){
        stop(sprintf(
            "No fragment overlaps the given regions in sample '%s'. Consider passing more regions or a larger 'windowSize'.",
            x ), call.=FALSE)
    }
    templateDf <- data.frame(
      relPos=seq_len(windowSize),
      fragmentSize=rep(seq_len(max(mcols(fragmentMidpoint)$fragmentSize)), each=windowSize))
    vMatrix <- as.data.frame(mcols(fragmentMidpoint))
    vMatrix <- aggregate(. ~ vMatrix$fragmentSize + vMatrix$relPos, data = vMatrix, FUN = length)
    vMatrix$fragmentSize <- NULL
    colnames(vMatrix)[colnames(vMatrix) == "vMatrix$fragmentSize"]<- "fragmentSize"
    colnames(vMatrix)[colnames(vMatrix) == "vMatrix$relPos"]<- "relPos"
    colnames(vMatrix)[ncol(vMatrix)] <- "n"
    vMatrix <- vMatrix[order(vMatrix$fragmentSize, vMatrix$relPos), ]
    rownames(vMatrix) <- NULL
    vMatrix <- merge(templateDf, vMatrix, by = c("fragmentSize", "relPos"), all.x = TRUE)

    vMatrix$n[is.na(vMatrix$n)] <- 0
    if( "Number of Reads" %in% colnames(colData(object)) ){
      vMatrix$totReads <- colData(object)[x,"Number of Reads"]
    }else{
      vMatrix$totReads <- sum( vMatrix$n )
    }
    if( is.null(sampleLabelColumn) ){
      vMatrix$samp <- x
    }else{
      vMatrix$samp <- colData(object)[x,sampleLabelColumn]
    }
    vMatrix
  }, BPPARAM=BPPARAM)
  vMatrices <- do.call(rbind, vMatrices)
  vMatrices$rpm <- (10^6) * vMatrices$n/vMatrices$totReads
  vMatrices$relPos <- vMatrices$relPos - (windowSize/2) + 1
  if( returnMatrix ){
    return(vMatrices)
  }else{
    if(rowScaled){
      vMatrices <- transform(vMatrices,
                             rpm = ave(rpm, samp, fragmentSize,
                             FUN = function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)))

      thres <- quantile(abs(vMatrices$rpm), zlim, na.rm=TRUE)
      vMatrices$rpm <-
        pmax(pmin(vMatrices$rpm, thres), -thres)
    }else{
      vMatrices$rpm <- pmin(vMatrices$rpm, quantile(vMatrices$rpm, zlim))
    }
    ## The `rpm` column holds z-scores when `rowScaled=TRUE`, so the legend
    ## has to follow the argument instead of always naming it RPM.
    fillName <- if( rowScaled ){
      "Row z-score"
    }else{
      "Fragment midpoints\n(RPM)"
    }
    scl <- scale_fill_gradient(
      low="#f7fbff", high="#08306b",
      na.value="#f7fbff", name=fillName)
    vMatrices %>%
      ggplot( aes( relPos, fragmentSize, fill=rpm ) ) +
      geom_tile() +
      coord_cartesian(ylim=c(0, 300)) +
      facet_wrap( ~samp ) +
      labs(x="Position relative to region centre (bp)",
           y="Fragment length (bp)",
           title="V-plot of fragment midpoints",
           caption=sprintf("Colour clipped at the %s%% quantile.",
                           format(100*zlim))) +
      scl +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
  }
}
