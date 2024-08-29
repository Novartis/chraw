#' Differential analyses of count experiments.
#'
#' @description This function inputs a ChrawExperiment object and the name of a count experiment represented
#' as a SummarizedExperiment within the ChrawExperiment object. Typically, users would to use the addCountExperiment
#' function to add these count experiments. These inputs are used for inference of differential
#' signals between conditions. Under the hood, DESeq2 is called and therefore some of the parameters are shared. The
#' outputs are stored as part of the rowRanges of the count experiment.
#'
#' @param object A ChrawExperiment object.
#' @param experimentName The name of the count experiment to test, normally added using the `addCountExperiment` function.
#' @param design A formula specifying the design to fit. See the examples or the ?DESeq documentation for more details.
#' @param contrasts A list specifying the contrasts to output. Each element of the list must be a character vector with three elements, c("name_of_the_variable", "test_condition", "baseline_condition"). See ?DESeq documentation for more details. If not specified, the last variable from the formula will be used and comparisons between all levels vs the first level will be returned.
#' @param replaceAll Logical indicating whether the existing results (if any) should be cleared and replace with new ones.
#' @param dropLevels Logical indicating whether unused levels should be kept in the testing step.
#' @param useColSizeFactors Logical. If `TRUE`, the numeric values from the column `sizeFactors` is used for library size normalization. Otherwise is calculated with the normal DESeq2 method.
#'
#' @return A ChrawExperiment object that include the results of the differential analysis.
#' @importFrom DESeq2 DESeq DESeqDataSetFromMatrix results
#' @importMethodsFrom S4Vectors elementMetadata
#' @export
testForDiffSignal <- function( object, experimentName, design, contrasts, replaceAll=FALSE, dropLevels=TRUE, useColSizeFactors=FALSE ){
    if( !is( object, "ChrawExperiment" ) ){
        stop( "The 'object' parameter must be a ChrawExperiment object" )
    }
    validObject( object )
    if( !experimentName %in% names(experiments(object)) ){
        stop(
            sprintf("The count experiment '%s' was not found in the ChrawExperiment. The value of the 'experimentName' parameter must match a count experiment name in the ChrawExperiment.",
                    experimentName))
    }
    if( !is( design, "formula") ){
        stop( "The 'design' parameter must be a formula." )
    }
    ## check whether the contrasts specify the right columns
    se <- experiments( object )[[experimentName]]
    if( !"counts" %in% assayNames( se ) ){
        stop( "No counts were found in the count experiment. Please make sure that the specified experiment has an assay names 'counts'" )
    }
    allVars <- all.vars( design )
    allVarsFlag <- allVars %in% colnames(colData(object))
    if( !all(allVarsFlag) ){
        stop(sprintf("All variables in the formula must be column names of colData(object). The following variable(s) from your formula are not in colnames(colData(object)):\n\t%s\n",
                     paste(allVars[!allVarsFlag], collapse="\n\t")))
    }
    ## check whether the contrasts specify the right columns
    if( !missing( contrasts ) ){
        if( !is( contrasts, "list" ) ){
            stop("The parameter 'contrasts' must be a named list.")
        }
        if( is.null(names(contrasts)) ){
            stop("The list passed to the parameter 'contrast' must be named.")
        }
        allContrastVars <- unique(unlist(lapply( contrasts, "[[", 1 )))
        allVarsFlag <- allContrastVars %in% allVars
        if( !all(allVarsFlag) ){
            stop(sprintf("The following variables specified in the contrasts parameters were not found in the design:\n\t%s\n",
                         paste(allContrastVars[!allVarsFlag], collapse="\n\t")))
        }
    }
    ## copy over the relevant columns from the colData of the ChrawExperiment to the count experiment
    for( i in allVars ){
        colData(se)[[i]] <- droplevels(factor(colData(object)[rownames(colData(se)),i]))
    }
    ## convert to DESeqDataSet for testing ##
    dsd <- DESeqDataSetFromMatrix(
        countData=as.matrix(assays(se)[["counts"]]),
        colData=colData(se)[,allVars,drop=FALSE],
        design=design )
    if( useColSizeFactors ){
        message("Using size factors defined in the ChrawExperiment object")
        dsd <- estimateSizeFactors(dsd)
        sizeFactors(dsd) <- colData(object)[rownames(colData(se)),"sizeFactors"]
    }
    if( missing( contrasts ) ){
        lvs <- levels(colData(dsd)[,allVars[length(allVars)]])
        lvs <- expand.grid(lvs[seq(2L, length(lvs))], lvs[1L])
        lvs <- as.matrix(lvs)
        lvs <- cbind( allVars[length(allVars)], lvs )
        lvs <- split( lvs, seq_len(nrow(lvs)) )
        names( lvs ) <- sprintf( "contrast%0.2d", seq_len( length( lvs ) ) )
        contrasts <- lvs
    }
    if( !dropLevels ){
        message("Fitting all data together...")
        dsd <- DESeq( dsd )
    }
    allResults <- lapply( names(contrasts), function(x){
        thisContrast <- contrasts[[x]]
        if( length(thisContrast) != 3L )
            stop(sprintf("Invalid contrast: %s, each contrast must be a character vector of length 3", paste(thisContrast, collapse=",")))
        if( !all(thisContrast[2:3] %in% colData(dsd)[[thisContrast[1L]]] ) )
            stop(sprintf("Invalid contrast %s, the reference and/or test levels were not found as values of the contrast variable",  paste(thisContrast, collapse=",")))
        if( dropLevels ){
            dsdSub <- dsd[,colData(dsd)[[thisContrast[1L]]] %in% thisContrast[2:3]]
            colData(dsdSub) <- droplevels(colData(dsdSub))
            message(sprintf("Fitting contrast %s...", paste(thisContrast, collapse=", ")))
            dsdSub <- DESeq(dsdSub)
        }else{
            dsdSub <- dsd
        }
        rs <- as( results( dsdSub, contrast=thisContrast ), "DataFrame" )
        elementMetadata( rs )$contrastName <- x
        elementMetadata( rs )$design <- paste(as.character(design(dsd)), collapse=" ")
        rs <- rs[,colnames(rs) %in% c("baseMean", "log2FoldChange", "pvalue", "padj")]
        colnames(rs) <- paste(elementMetadata( rs )$contrastName, colnames(rs), sep="_")
        rs
    } )
    allResults <- do.call( cbind, allResults )
    elementMetadata(allResults)$type <- "differential_test"
    ## appending the results to the current rowData
    ## replace or substitute results according to the parameters specified
    if( any(elementMetadata(rowData(se))$type == "differential_test") ){
        if( replaceAll ){
            rowData(se) <- rowData(se)[,!elementMetadata(rowData(se))$type == "differential_test"]
        }
        existsFlag <- elementMetadata(rowData(se))$contrastName %in% names(contrasts) &
                                                 elementMetadata(rowData(se))$type == "differential_test"
        if( any( existsFlag ) ){
            warning("The ChrawExperiment object has already results with the same contrast names, replacing these...")
            rowData(se) <- rowData(se)[,!existsFlag]
        }
    }
    rowData( se ) <- cbind( rowData(se), allResults )
    experiments(object)[[experimentName]] <- se
    object
}

#' Returns all the results of diffential analyses stored in a ChrawExperiment object.
#'
#' @description This function loops over the `experiments()` that were added using the
#' `addCountData()` function to the ChrawExperiment object, to look for results of
#' differential testing generated using the `testForDiffSignal()`. It returns a
#' `data.frame()` that summarizes the results stored in the ChrawExperiment object.
#'
#' @param object A ChrawExperiment object.
#' @return A `data.frame()` summarizing the ChrawExperiment object that include the results of the differential analysis.
#'
#' @export
getDiffResultSummary <- function( object ){
    if( !is(object, "ChrawExperiment") ){
        stop("The parameter 'object' must be a ChrawExperiment object")
    }
    lookExps <- vapply( experiments(object), function(x){is(x, "SummarizedExperiment")}, logical(1) )
    if( !any( lookExps ) ){
        stop("No SummarizedExperiments were found in the `experiments()` of the ChrawExperiment object. Check ?addCountExperiment to add these.")
    }
    allResults <- lapply( names(experiments(object)[lookExps]), function(x){
        xname <- x
        x <- experiments(object)[[x]]
        resultsColumns <- elementMetadata(rowData(x))$type == "differential_test"
        if( sum( resultsColumns, na.rm=TRUE ) > 0 ){
            resSummary <- elementMetadata(rowData(x)[,which(resultsColumns)])
            resSummary <- resSummary[grepl("log2 fold change", resSummary$description),]
            rs <- data.frame(
                experimentName=xname,
                contrastName=resSummary$contrastName,
                design=resSummary$design,
                description=gsub("log2 fold change \\(MLE\\): ", "\\1", resSummary$description) )
            rs$testVariable <- gsub("(\\S+)\\s.*", "\\1", rs$description)
            rs$testLevel <- gsub("\\S+\\s(\\S+) vs \\S+", "\\1", rs$description)
            rs$referenceLevel <- gsub("\\S+\\s\\S+ vs (\\S+)", "\\1", rs$description)
            rs$description <- NULL
            return(rs)
        }else{
            return(NULL)
        }
    })
    allResults <- do.call(rbind, allResults)
    if( is.null( allResults ) ){
        stop("No results were found in the ChrawExperiment object, nothing to return")
    }
    allResults
}

#' Extracts the results of the differential analysis from a ChrawExperiment object
#'
#' @description Having a ChrawExperiment object where results from a differential
#' analysis were run and stored using the `testForDiffSignal()` object, these function
#' extract the results of the differential analysis specified in the parameters described
#' in this man page. For a list of the stored results, consider running the `getDiffResultSummary()`
#' function.
#'
#' @param object A ChrawExperiment object.
#' @param experimentName The name of an count experiment added using the `addCountExperiment()` function.
#' @param contrastName A name of the contrast name to extract. For a list of the contrast names, consider running `getDiffResultSummary()`. If missing, all available results will be returned.
#' @param outputFormat A character string with one of the following values: "GRanges" (for a GRanges object), "df_wide" (wide data frame), or "df_long" (long data frame).
#' @return An object containing the results of the differential analysis. The specific format is specified in the `outputFormat=` parameter.
#'
#' @export
pullDiffResults <- function( object, experimentName, contrastName=NULL, outputFormat="GRanges" ){
    if( !is( object, "ChrawExperiment" ) ){
        stop("The parameter 'object' must be a ChrawExperiment object")
    }
    if( !experimentName %in% names( experiments( object ) ) ){
        stop(sprintf("The experiment '%s' was not found in the ChrawExperiment object. Check `names(experiments(object))` for valid values.", experimentName))
    }
    if( !outputFormat %in% c("GRanges", "df_wide", "df_long") ){
        stop("The outputFormat values must be either 'GRanges', 'df_wide' or 'df_long'.")
    }
    se <- experiments( object )[[experimentName]]
    if( !is.null(contrastName) ){
        if( !all(contrastName %in% elementMetadata( rowData( se ) )$contrastName) ){
            stop(sprintf(
                "The contrast(s) '%s' was/were not found in the ChrawExperiment object. Please run `getDiffResultSummary()` to see the available contrast names.",
                paste( contrastName[!contrastName %in% elementMetadata( rowData( se ) )$contrastName], collapse=", " ) ) )
        }
    }else{
        ## if missing contrastName parameter, specify all contrasts
        contrastName <- unique(elementMetadata( rowData( se ) )$contrastName)
    }
    ## extracts the desired outputs
    diffInfo <- rowData(se)[,elementMetadata( rowData( se ) )$contrastName %in% contrastName]
    if( ncol(diffInfo) == 0 ){
        stop("No differential results were found. Consider running `testForDiffSignal()` to add differential results.")
    }
    rr <- rowRanges(se)
    mcols(rr) <- NULL
    mcols(rr) <- diffInfo
    ## reformat output according to `outputFormat=` parameter.
    if( outputFormat == "GRanges" ){
        return(rr)
    }else if( outputFormat == "df_wide" ){
        dfRes <- cbind(
            data.frame(
                regionID=names(rr),
                genomicCoordinates=as.character(rr) ),
            as.data.frame( mcols(rr) ) )
        return(dfRes)
    }else if( outputFormat == "df_long" ){
        spCols <- split( seq_len(ncol(mcols(rr))), elementMetadata( mcols(rr) )$contrastName )
        dfRes <- lapply( names(spCols), function(x){
            resSub <- mcols(rr)[,spCols[[x]]]
            colnames(resSub) <- gsub(paste0(x, "_"), "", colnames(resSub))
            as.data.frame(cbind(
                contrastName=x,
                regionID=names(rr),
                genomicCoordinates=as.character(rr), resSub ) )
        } )
        dfRes <- do.call(rbind, dfRes)
        return(dfRes)
    }
}
