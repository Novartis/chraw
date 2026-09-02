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
#' This function defines the S4 method indexBam for the object
#'   ChrawExperiment.
#' This method locates the bam file of a ChrawExperiment and creates the
#'   indexes
#' for the bam files that are not indexed.
#'
#' @docType methods
#' @name indexBam
#' @rdname indexBam
#' @aliases indexBam indexBam,ChrawExperiment-method
#'
#' @param files A ChrawExperiment object.
#' @param dedup Logical indicating whether to index the deduplicated bam
#'   files or the full bam files.
#' @param BPPARAM A BiocParallel instance.
#' @param force Logical indicating whether indexes should be re-created in
#'   case they exist.
#'
#' @import BiocParallel
#' @import Rsamtools
#' @return Invisibly returns `NULL`. Called for the side effect of
#' writing a `.bai` index next to each bam file that lacks one.
#'
#' @examples
#' data(ce_chipseq)
#' ce_chipseq <- rewrite_paths(ce_chipseq)
#' 
#' indexBam( ce_chipseq )
#'
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
#' @aliases getGenomicBins getGenomicBins,BSgenome-method
#'   getGenomicBins,ChrawExperiment-method
#'
#' @param object A BSgenome object or a ChrawExperiment object
#' @param binSize A numeric value specifying the size of the bins desired.
#' @param onlyStandardChromosomes Logical indicating whether to keep only
#'   standard chromosomes. If 'FALSE' chromosome patches and haplotypes are
#'   dropped, as well as the mitochrondrial and 'Y' chromosomes.
#'
#' @return A [GenomicRanges::GRanges] object with the genomic bins.
#'
#' @examples
#' data(ce_examples)
#' ce_examples <- rewrite_paths(ce_examples)
#' 
#' bins <- getGenomicBins( ce_examples, binSize = 10^6 )
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
#' @examples
#' data(ce_chipseq)
#' ce_chipseq <- rewrite_paths(ce_chipseq)
#' indexBam( ce_chipseq )
#' 
#' fragLengths <- computeFragLengthDist(
#'     colData(ce_chipseq)$bamFile[1],
#'     param = csaw::readParam(pe = "both", restrict = "chr6") )
#' head( fragLengths )
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
        if( length(args) == 1 ){
            input <- args[[1L]]
            if( !is(input, "ChrawExperiment")){
                return(callNextMethod(x, ..., sampleMap=sampleMap, mapFrom=mapFrom ))
            }

            refGenome <- unique(c(referenceGenome(x), referenceGenome(input)))
            if( length( refGenome ) > 1 ){
                stop("Merging 'ChrawExperiment' objects from different species is not supported")
            }
            pipe <- unique(c(pipeline(x), pipeline(input)))
            if( length( pipe ) > 1 ){
                warning("Merging 'ChrawExperiment' objects from different pipelines, QC functionality will be disabled")
                pipe <- "public"
            }

            existsMatrix <- cbind(experiments(x)[["ExistingSampleFlag"]],
                                  experiments(input)[["ExistingSampleFlag"]])
            tmpName1 <- as.vector(randomStrings(1, len=10))
            tmpName2 <- as.vector(randomStrings(1, len=10))
            x <- addExperiment( x, experiments(x)[["ExistingSampleFlag"]], tmpName1 )
            input <- addExperiment( input, experiments(input)[["ExistingSampleFlag"]], tmpName2 )
            suppressMessages(experiments(x)[["ExistingSampleFlag"]] <- NULL)
            suppressMessages(experiments(input)[["ExistingSampleFlag"]] <- NULL)
            concatCe <- callNextMethod(x, input)
            concatCe <- addExperiment( concatCe, existsMatrix, "ExistingSampleFlag" )
            suppressMessages(experiments(concatCe)[[tmpName1]] <- NULL)
            suppressMessages(experiments(concatCe)[[tmpName2]] <- NULL)
            cr <- new("ChrawExperiment",
                      concatCe,
                      referenceGenome=refGenome,
                      pipeline=pipe )
            return(cr)
        }else if( length(args) > 1 ){
            for( i in seq_along(args) ){
                x <- c( x, args[[i]], sampleMap=sampleMap, mapFrom=mapFrom )
            }
            return(x)
        }
    })

#' Accessors for the metadata of a ChrawExperiment object
#'
#' @description `referenceGenome()` returns the reference genome assembly that
#' the experiments of a [ChrawExperiment-class] object were aligned to, and
#' `pipeline()` returns the name of the pipeline that produced them. Both are
#' recorded when the object is built and are used across the package to select
#' the matching annotation resources.
#'
#' @docType methods
#' @name referenceGenome
#' @rdname referenceGenome
#' @aliases referenceGenome referenceGenome,ChrawExperiment-method
#' pipeline pipeline,ChrawExperiment-method
#'
#' @param x A ChrawExperiment object.
#' @param object A ChrawExperiment object.
#'
#' @return A character vector of length 1 with the reference genome assembly
#' (`referenceGenome()`) or the pipeline name (`pipeline()`).
#'
#' @examples
#' data(ce_examples)
#'
#' referenceGenome( ce_examples )
#' pipeline( ce_examples )
#'
#' @importFrom BSgenome referenceGenome
#' @export
setMethod( "referenceGenome", signature( x="ChrawExperiment" ),
          function( x ) x@referenceGenome )

#' @rdname referenceGenome
#' @export
setMethod( "pipeline", signature( object="ChrawExperiment" ),
          function( object ) object@pipeline )
