indexCromwellBams <- function( files, dedup=FALSE, BPPARAM=SerialParam(), force=FALSE ){
    if( dedup )
        colData(files)$bamFile <- getDedupAlignments( colData(files)$bamFile )
    bfls <- BamFileList(colData( files )$bamFile)
    if( !force )
        bfls <- bfls[lengths(lapply(bfls, "[[", "index")) == 0]
    if( length( bfls ) > 0 ){
        indexes <- bplapply(
            bfls,
            function(x){
                indexBam( x$path )
            }, BPPARAM=BPPARAM )
        logicalIdx <- file.exists(unlist(indexes))
        if( !all(logicalIdx) ){
            warning(sprintf("Index files could not be created for the following bam files:\n\t%s",
                            paste(names(bfls)[!logicalIdx], collapse="\n\t")))
        }
    }
    invisible()
}

#' Index the bam files from a ChrawExperiment object
#'
#' This function defines the S4 method indexBam for the object ChrawExperiment.
#' This method locates the bam file of a ChrawExperiment and creates the indexes
#' for the bam files that are not indexed.
#'
#' @docType methods
#' @name indexBam
#' @rdname indexBam
#' @aliases indexBam indexBam,ChrawExperiment-method
#'
#' @param files A ChrawExperiment object.
#' @param dedup Logical indicating whether to index the deduplicated bam files or the full bam files.
#' @param BPPARAM A BiocParallel instance.
#' @param force Logical indicating whether indexes should be re-created in case they exist.
#'
#' @import BiocParallel
#' @import Rsamtools
#' @importMethodsFrom Rsamtools indexBam
#' @export
setMethod( "indexBam", signature(files="ChrawExperiment"), indexCromwellBams )

getGenomicBins.bsgenome <- function( object, binSize, onlyStandardChromosomes=TRUE ){
    if(!is.numeric(binSize))
        stop("The parameter 'binSize' must be a numeric value")
    seqLengths <- seqlengths( object )
    compartmentGR <- do.call(rbind, lapply( names(seqLengths), function(chr){
        st <- seq( 1, seqLengths[chr], binSize )
        end <- st + binSize - 1
        end <- pmin(end, seqLengths[chr])
        data.frame(chr=chr, start=st, end=end)
    } ) )
    compartmentGR <- makeGRangesFromDataFrame(compartmentGR)
    mcols(compartmentGR) <- NULL
    if( onlyStandardChromosomes ){
        compartmentGR <- keepStandardChromosomes( compartmentGR, pruning.mode="coarse" )
        compartmentGR <- dropSeqlevels( compartmentGR, c("chrM", "chrY"), pruning.mode="coarse")
    }
    compartmentGR
}

getGenomicBins.ChrawExperiment <- function( object, binSize, onlyStandardChromosomes=TRUE ){
    validObject(object)
    getGenomicBins( selectBSgenome( object ), binSize, onlyStandardChromosomes )
}

#' Partition a reference genome in bins
#'
#' This method inputs either a ChrawExperiment or a BSgenome and returns
#' a GenomicRanges object with equally sized genomic bins across the genome.
#'
#' @docType methods
#' @name getGenomicBins
#' @rdname getGenomicBins
#' @aliases getGenomicBins getGenomicBins,BSgenome-method getGenomicBins,ChrawExperiment-method
#'
#' @param object A BSgenome object or a ChrawExperiment object
#' @param binSize A numeric value specifying the size of the bins desired.
#' @param onlyStandardChromosomes Logical indicating whether to keep only standard chromosomes. If 'FALSE' chromosome patches and haplotypes are dropped, as well as the mitochrondrial and 'Y' chromosomes.
#'
#' @import GenomicRanges
#' @export
setMethod( getGenomicBins, signature( object = "BSgenome" ), getGenomicBins.bsgenome )

#' @name getGenomicBins
#' @rdname getGenomicBins
#' @export
setMethod( getGenomicBins, signature( object = "ChrawExperiment" ), getGenomicBins.ChrawExperiment )


#' Compute fragment length distribution from a paired-end *.bam file.
#'
#' @param bam A character string specifying the file paths to bam file.
#' @param ... Additional parameters passed to `getPESizes()`.
#'
#' @importFrom csaw getPESizes
#'
#' @return Data frame with fragment lenghts and thei frequency
#' @export
#'
computeFragLengthDist <- function( bam, ... ) {
    frag.lengths <- getPESizes( bam, ... )
    frag.lengths <- table(frag.lengths$sizes)
    frag.lengths <- data.frame(bamReads = bam,
                               fraglen = as.numeric(names(frag.lengths)),
                               freq = as.numeric(frag.lengths)/sum(as.numeric(frag.lengths)))
    return(frag.lengths)
}

#' @rdname addCountExperiment
#'
#' @param x A `ChrawExperiment` object.
#' @param sampleMap \code{c} method: a \code{sampleMap} \code{list} or
#' \code{DataFrame} to guide merge
#' @param mapFrom Either a \code{logical}, \code{character}, or \code{integer}
#' vector indicating the experiment(s) that have an identical colname order as
#' the experiment input(s). If using a character input, the name must match
#' exactly.
#' @export
setMethod(
    "c", c(x="ChrawExperiment"),
    function( x, ..., sampleMap = NULL, mapFrom = NULL) {
        args <- list(...)
        ## print(length(args))
        if( length(args) == 1 ){
            ## print("entro")
            input <- args[[1L]]
            ## print(class(input))
            if( !is(input, "ChrawExperiment")){
                return(callNextMethod(x, ..., sampleMap=sampleMap, mapFrom=mapFrom ))
            }

            referenceGenome <- unique(c(x@referenceGenome, input@referenceGenome))
            if( length( referenceGenome ) > 1 ){
                stop("Merging 'ChrawExperiment' objects from different species is not supported")
            }
            pipeline <- unique(c(x@pipeline, input@pipeline))
            if( length( pipeline ) > 1 ){
                warning("Merging 'ChrawExperiment' objects from different pipelines, QC functionality will be disabled")
                pipeline <- "public"
            }


            existsMatrix <- cbind(experiments(x)[["ExistingSampleFlag"]], experiments(input)[["ExistingSampleFlag"]])
            expNames1 <- names(experiments(x))
            expNames2 <- names(experiments(input))
            tmpName1 <- as.vector(randomStrings(1, len=10))
            tmpName2 <- as.vector(randomStrings(1, len=10))
            #experiments(x)[[tmpName1]] <- experiments(x)[["ExistingSampleFlag"]]
            #experiments(input)[[tmpName2]] <- experiments(input)[["ExistingSampleFlag"]]
            x <- addExperiment( x, experiments(x)[["ExistingSampleFlag"]], tmpName1 )
            input <- addExperiment( input, experiments(input)[["ExistingSampleFlag"]], tmpName2 )
            suppressMessages(experiments(x)[["ExistingSampleFlag"]] <- NULL)
            suppressMessages(experiments(input)[["ExistingSampleFlag"]] <- NULL)
            ##concatCe <- c(as(x, "MultiAssayExperiment", as(input, "MultiAssayExperiment")))
            concatCe <- callNextMethod(x, input)
            concatCe <- addExperiment( concatCe, existsMatrix, "ExistingSampleFlag" )
            suppressMessages(experiments(concatCe)[[tmpName1]] <- NULL)
            suppressMessages(experiments(concatCe)[[tmpName2]] <- NULL)
            cr <- new("ChrawExperiment",
                      concatCe,
                      referenceGenome=referenceGenome,
                      pipeline=pipeline )
            return(cr)
        }else if( length(args) > 1 ){
            for( i in seq_len(length(args)) ){
                x <- c( x, args[[i]], sampleMap=sampleMap, mapFrom=mapFrom )
            }
            return(x)
        }
    })
